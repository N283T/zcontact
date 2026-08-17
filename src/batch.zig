//! Deterministic, file-parallel batch processing for flat structure directories.

const std = @import("std");
const contact = @import("contact.zig");
const gzip = @import("gzip.zig");
const output = @import("output.zig");
const parser = @import("parser.zig");

pub const Config = struct {
    input_dir: []const u8,
    output_dir: []const u8,
    manifest_path: ?[]const u8 = null,
    mode: contact.Mode = .residue,
    cutoff: f64 = 4.0,
    select1: []const u8 = "polymer,heavy",
    select2: []const u8 = "polymer,heavy",
    scope: contact.Scope = .inter_residue,
    model: u32 = 1,
    threads: u32 = 0,
    skip_existing: bool = false,
    overwrite: bool = false,
    quiet: bool = false,
    tool_version: []const u8 = "unknown",
};

pub const Status = enum { ok, skipped, err };

pub const FileResult = struct {
    status: Status = .err,
    atom_count: usize = 0,
    residue_count: usize = 0,
    contact_count: usize = 0,
    input_bytes: u64 = 0,
    output_bytes: u64 = 0,
    elapsed_ns: u64 = 0,
    error_name: []const u8 = "",
};

pub const Summary = struct {
    total: usize,
    successful: usize,
    skipped: usize,
    failed: usize,
    contacts: usize,
    wall_time_ns: u64,
};

fn isStructureFile(name: []const u8) bool {
    return std.mem.endsWith(u8, name, ".pdb") or std.mem.endsWith(u8, name, ".pdb.gz") or
        std.mem.endsWith(u8, name, ".ent") or std.mem.endsWith(u8, name, ".ent.gz") or
        std.mem.endsWith(u8, name, ".cif") or std.mem.endsWith(u8, name, ".cif.gz") or
        std.mem.endsWith(u8, name, ".mmcif") or std.mem.endsWith(u8, name, ".mmcif.gz");
}

pub fn scanDirectory(allocator: std.mem.Allocator, path: []const u8) ![][]const u8 {
    const io = std.Io.Threaded.global_single_threaded.io();
    var dir = try std.Io.Dir.cwd().openDir(io, path, .{ .iterate = true });
    defer dir.close(io);
    var names = std.ArrayListUnmanaged([]const u8).empty;
    errdefer {
        for (names.items) |name| allocator.free(name);
        names.deinit(allocator);
    }
    var iterator = dir.iterate();
    while (try iterator.next(io)) |entry| {
        if (entry.kind != .file or entry.name.len == 0 or entry.name[0] == '.' or !isStructureFile(entry.name)) continue;
        if (std.mem.indexOfAny(u8, entry.name, "\t\r\n") != null) return error.UnsafeFilename;
        try names.append(allocator, try allocator.dupe(u8, entry.name));
    }
    const result = try names.toOwnedSlice(allocator);
    std.mem.sort([]const u8, result, {}, struct {
        fn lessThan(_: void, a: []const u8, b: []const u8) bool {
            return std.mem.lessThan(u8, a, b);
        }
    }.lessThan);
    return result;
}

pub fn freeFileList(allocator: std.mem.Allocator, files: [][]const u8) void {
    for (files) |name| allocator.free(name);
    allocator.free(files);
}

fn outputName(allocator: std.mem.Allocator, filename: []const u8) ![]u8 {
    return std.fmt.allocPrint(allocator, "{s}.contacts.tsv", .{std.fs.path.basename(filename)});
}

fn sha256Bytes(data: []const u8) [64]u8 {
    var digest: [32]u8 = undefined;
    std.crypto.hash.sha2.Sha256.hash(data, &digest, .{});
    return std.fmt.bytesToHex(digest, .lower);
}

fn sha256File(io: std.Io, path: []const u8) ![64]u8 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    var reader_buffer: [64 * 1024]u8 = undefined;
    var read_buffer: [64 * 1024]u8 = undefined;
    var reader = file.reader(io, &reader_buffer);
    var hasher = std.crypto.hash.sha2.Sha256.init(.{});
    while (true) {
        const count = try reader.interface.readSliceShort(&read_buffer);
        if (count == 0) break;
        hasher.update(read_buffer[0..count]);
    }
    var digest: [32]u8 = undefined;
    hasher.final(&digest);
    return std.fmt.bytesToHex(digest, .lower);
}

fn fileSize(io: std.Io, path: []const u8) !u64 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    return (try file.stat(io)).size;
}

fn tempPath(allocator: std.mem.Allocator, io: std.Io, path: []const u8) ![]u8 {
    return std.fmt.allocPrint(allocator, "{s}.tmp-{d}-{d}", .{ path, std.Thread.getCurrentId(), std.Io.Timestamp.now(io, .awake).nanoseconds });
}

fn markerPath(allocator: std.mem.Allocator, output_path: []const u8) ![]u8 {
    return std.fmt.allocPrint(allocator, "{s}.meta.tsv", .{output_path});
}

fn canonicalPath(allocator: std.mem.Allocator, io: std.Io, path: []const u8) ![]u8 {
    const real = std.Io.Dir.cwd().realPathFileAlloc(io, path, allocator) catch |err| switch (err) {
        error.FileNotFound => return std.fs.path.resolve(allocator, &.{path}),
        else => return err,
    };
    defer allocator.free(real);
    return allocator.dupe(u8, real);
}

fn markerMatches(
    allocator: std.mem.Allocator,
    config: *const Config,
    marker_path: []const u8,
    input_sha: [64]u8,
    output_sha: [64]u8,
    input_bytes: u64,
    output_bytes: u64,
) !?FileResult {
    const io = std.Io.Threaded.global_single_threaded.io();
    const data = std.Io.Dir.cwd().readFileAlloc(io, marker_path, allocator, .limited(16 * 1024)) catch return null;
    defer allocator.free(data);
    var lines = std.mem.splitScalar(u8, data, '\n');
    _ = lines.next() orelse return null;
    const row = lines.next() orelse return null;
    var fields: [14][]const u8 = undefined;
    var count: usize = 0;
    var parts = std.mem.splitScalar(u8, row, '\t');
    while (parts.next()) |part| {
        if (count >= fields.len) return null;
        fields[count] = part;
        count += 1;
    }
    if (count != fields.len or !std.mem.eql(u8, fields[0], "1") or
        !std.mem.eql(u8, fields[1], @tagName(config.mode)) or
        !std.mem.eql(u8, fields[2], config.select1) or !std.mem.eql(u8, fields[3], config.select2) or
        !std.mem.eql(u8, fields[4], @tagName(config.scope)) or
        (std.fmt.parseInt(u32, fields[5], 10) catch return null) != config.model or
        (std.fmt.parseFloat(f64, fields[6]) catch return null) != config.cutoff or
        !std.mem.eql(u8, fields[7], &input_sha) or !std.mem.eql(u8, fields[8], &output_sha) or
        (std.fmt.parseInt(u64, fields[12], 10) catch return null) != input_bytes or
        (std.fmt.parseInt(u64, fields[13], 10) catch return null) != output_bytes) return null;
    return .{
        .status = .skipped,
        .atom_count = std.fmt.parseInt(usize, fields[9], 10) catch return null,
        .residue_count = std.fmt.parseInt(usize, fields[10], 10) catch return null,
        .contact_count = std.fmt.parseInt(usize, fields[11], 10) catch return null,
        .input_bytes = input_bytes,
        .output_bytes = output_bytes,
    };
}

fn writeMarker(allocator: std.mem.Allocator, config: *const Config, path: []const u8, result: FileResult, input_sha: [64]u8, output_sha: [64]u8) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const temp_path = try tempPath(allocator, io, path);
    defer allocator.free(temp_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, temp_path) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, temp_path, .{ .exclusive = true });
        defer file.close(io);
        var buffer: [4096]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try writer.interface.writeAll("schema\tmode\tselect1\tselect2\tscope\tmodel\tcutoff_A\tinput_sha256\toutput_sha256\tatoms\tresidues\tcontacts\tinput_bytes\toutput_bytes\n");
        try writer.interface.print("1\t{s}\t{s}\t{s}\t{s}\t{d}\t{d}\t{s}\t{s}\t{d}\t{d}\t{d}\t{d}\t{d}\n", .{ @tagName(config.mode), config.select1, config.select2, @tagName(config.scope), config.model, config.cutoff, &input_sha, &output_sha, result.atom_count, result.residue_count, result.contact_count, result.input_bytes, result.output_bytes });
        try writer.interface.flush();
    }
    try std.Io.Dir.renameAbsolute(temp_path, path, io);
}

fn processOne(allocator: std.mem.Allocator, config: *const Config, filename: []const u8) !FileResult {
    const io = std.Io.Threaded.global_single_threaded.io();
    const input_path = try std.fs.path.join(allocator, &.{ config.input_dir, filename });
    defer allocator.free(input_path);
    const out_name = try outputName(allocator, filename);
    defer allocator.free(out_name);
    const output_path = try std.fs.path.join(allocator, &.{ config.output_dir, out_name });
    defer allocator.free(output_path);
    const marker_path = try markerPath(allocator, output_path);
    defer allocator.free(marker_path);

    const exists = blk: {
        std.Io.Dir.cwd().access(io, output_path, .{}) catch |err| switch (err) {
            error.FileNotFound => break :blk false,
            else => return err,
        };
        break :blk true;
    };
    if (exists) {
        if (!config.skip_existing and !config.overwrite) return error.OutputExists;
    }

    const started = std.Io.Timestamp.now(io, .awake);
    const input_bytes = try fileSize(io, input_path);
    const source = try gzip.readFile(allocator, input_path, gzip.default_max_size);
    defer allocator.free(source);
    const input_sha = sha256Bytes(source);
    if (exists and config.skip_existing) {
        const output_bytes = try fileSize(io, output_path);
        const output_sha = try sha256File(io, output_path);
        if (try markerMatches(allocator, config, marker_path, input_sha, output_sha, input_bytes, output_bytes)) |matched| {
            var valid = matched;
            valid.elapsed_ns = @intCast(started.untilNow(io, .awake).nanoseconds);
            return valid;
        }
    }
    const format = try parser.detectFormat(input_path, source);
    var structure = try parser.parse(allocator, source, format, config.model);
    defer structure.deinit(allocator);
    var contacts = try contact.calculate(allocator, &structure, config.mode, config.cutoff, config.select1, config.select2, config.scope);
    defer contacts.deinit(allocator);

    const temp_path = try tempPath(allocator, io, output_path);
    defer allocator.free(temp_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, temp_path) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, temp_path, .{ .exclusive = true });
        defer file.close(io);
        var buffer: [64 * 1024]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try output.writeResult(&writer.interface, &structure, &contacts);
        try writer.interface.flush();
    }
    try std.Io.Dir.renameAbsolute(temp_path, output_path, io);
    const output_bytes = try fileSize(io, output_path);
    var residue_count: usize = 0;
    for (structure.atoms.items) |atom| residue_count = @max(residue_count, @as(usize, atom.residue_index) + 1);
    const output_sha = try sha256File(io, output_path);
    const result = FileResult{
        .status = .ok,
        .atom_count = structure.atoms.items.len,
        .residue_count = residue_count,
        .contact_count = output.count(&contacts),
        .input_bytes = input_bytes,
        .output_bytes = output_bytes,
        .elapsed_ns = @intCast(started.untilNow(io, .awake).nanoseconds),
    };
    try writeMarker(allocator, config, marker_path, result, input_sha, output_sha);
    return result;
}

const WorkerContext = struct {
    config: *const Config,
    files: []const []const u8,
    results: []FileResult,
    next: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    completed: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
};

fn worker(context: *WorkerContext) void {
    const io = std.Io.Threaded.global_single_threaded.io();
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();
    while (true) {
        const index = context.next.fetchAdd(1, .monotonic);
        if (index >= context.files.len) return;
        // Do not retain the peak allocation of a very large structure for the
        // lifetime of the worker. Batch parallelism already multiplies the
        // whole-file parser's memory use.
        _ = arena.reset(.free_all);
        const started = std.Io.Timestamp.now(io, .awake);
        context.results[index] = processOne(arena.allocator(), context.config, context.files[index]) catch |err| .{
            .status = .err,
            .input_bytes = blk: {
                const path = std.fs.path.join(arena.allocator(), &.{ context.config.input_dir, context.files[index] }) catch break :blk 0;
                break :blk fileSize(io, path) catch 0;
            },
            .elapsed_ns = @intCast(started.untilNow(io, .awake).nanoseconds),
            .error_name = @errorName(err),
        };
        const done = context.completed.fetchAdd(1, .monotonic) + 1;
        if (!context.config.quiet and (done % 100 == 0 or done == context.files.len))
            std.debug.print("\rProcessed {d}/{d}", .{ done, context.files.len });
    }
}

fn writeManifest(allocator: std.mem.Allocator, config: *const Config, files: []const []const u8, results: []const FileResult) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const owned_path = if (config.manifest_path) |path|
        try allocator.dupe(u8, path)
    else
        try std.fs.path.join(allocator, &.{ config.output_dir, "manifest.tsv" });
    defer allocator.free(owned_path);
    const temp_path = try tempPath(allocator, io, owned_path);
    defer allocator.free(temp_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, temp_path) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, temp_path, .{ .exclusive = true });
        defer file.close(io);
        var buffer: [64 * 1024]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try writer.interface.writeAll("input\toutput\tstatus\tatoms\tresidues\tcontacts\tinput_bytes\toutput_bytes\telapsed_ms\terror\n");
        for (files, results) |name, result| {
            const out_name = try outputName(allocator, name);
            defer allocator.free(out_name);
            try writer.interface.print("{s}\t{s}\t{s}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d:.3}\t{s}\n", .{
                name,
                out_name,
                @tagName(result.status),
                result.atom_count,
                result.residue_count,
                result.contact_count,
                result.input_bytes,
                result.output_bytes,
                @as(f64, @floatFromInt(result.elapsed_ns)) / 1_000_000.0,
                result.error_name,
            });
        }
        try writer.interface.flush();
    }
    try std.Io.Dir.renameAbsolute(temp_path, owned_path, io);
}

fn prepareRunMetadata(allocator: std.mem.Allocator, config: *const Config) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const path = try std.fs.path.join(allocator, &.{ config.output_dir, "run.tsv" });
    defer allocator.free(path);
    var expected_writer: std.Io.Writer.Allocating = .init(allocator);
    defer expected_writer.deinit();
    try expected_writer.writer.print(
        "key\tvalue\nschema_version\t1\nzcontact_version\t{s}\ninput_dir\t{s}\nmode\t{s}\ncutoff_A\t{d}\nselect1\t{s}\nselect2\t{s}\nscope\t{s}\nmodel\t{d}\n",
        .{ config.tool_version, config.input_dir, @tagName(config.mode), config.cutoff, config.select1, config.select2, @tagName(config.scope), config.model },
    );
    const expected = expected_writer.writer.buffered();
    const existing = std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(64 * 1024)) catch |err| switch (err) {
        error.FileNotFound => null,
        else => return err,
    };
    defer if (existing) |data| allocator.free(data);
    if (config.skip_existing) {
        const actual = existing orelse return error.ResumeMetadataMissing;
        if (!std.mem.eql(u8, actual, expected)) return error.ResumeConfigurationMismatch;
        return;
    }
    if (existing != null and !config.overwrite) return error.BatchDirectoryInitialized;

    const temp_path = try tempPath(allocator, io, path);
    defer allocator.free(temp_path);
    errdefer std.Io.Dir.cwd().deleteFile(io, temp_path) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, temp_path, .{ .exclusive = true });
        defer file.close(io);
        var buffer: [4096]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try writer.interface.writeAll(expected);
        try writer.interface.flush();
    }
    try std.Io.Dir.renameAbsolute(temp_path, path, io);
}

fn validateManifestPath(allocator: std.mem.Allocator, config: *const Config, files: []const []const u8) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const default_manifest = try std.fs.path.join(allocator, &.{ config.output_dir, "manifest.tsv" });
    defer allocator.free(default_manifest);
    const manifest = try canonicalPath(allocator, io, config.manifest_path orelse default_manifest);
    defer allocator.free(manifest);
    const run_path_raw = try std.fs.path.join(allocator, &.{ config.output_dir, "run.tsv" });
    defer allocator.free(run_path_raw);
    const run_path = try canonicalPath(allocator, io, run_path_raw);
    defer allocator.free(run_path);
    if (std.mem.eql(u8, manifest, run_path)) return error.ManifestPathCollision;
    for (files) |name| {
        const out_name = try outputName(allocator, name);
        defer allocator.free(out_name);
        const output_path_raw = try std.fs.path.join(allocator, &.{ config.output_dir, out_name });
        defer allocator.free(output_path_raw);
        const output_path = try canonicalPath(allocator, io, output_path_raw);
        defer allocator.free(output_path);
        if (std.mem.eql(u8, manifest, output_path)) return error.ManifestPathCollision;
        const output_marker_raw = try markerPath(allocator, output_path);
        defer allocator.free(output_marker_raw);
        const output_marker = try canonicalPath(allocator, io, output_marker_raw);
        defer allocator.free(output_marker);
        if (std.mem.eql(u8, manifest, output_marker)) return error.ManifestPathCollision;
        const input_path_raw = try std.fs.path.join(allocator, &.{ config.input_dir, name });
        defer allocator.free(input_path_raw);
        const input_path = try canonicalPath(allocator, io, input_path_raw);
        defer allocator.free(input_path);
        if (std.mem.eql(u8, manifest, input_path)) return error.ManifestPathCollision;
    }
    if (config.manifest_path != null and !config.skip_existing and !config.overwrite)
        std.Io.Dir.cwd().access(io, manifest, .{}) catch |err| switch (err) {
            error.FileNotFound => return,
            else => return err,
        }
    else
        return;
    return error.ManifestExists;
}

pub fn run(allocator: std.mem.Allocator, config: *const Config) !Summary {
    const io = std.Io.Threaded.global_single_threaded.io();
    const started = std.Io.Timestamp.now(io, .awake);
    const files = try scanDirectory(allocator, config.input_dir);
    defer freeFileList(allocator, files);
    if (files.len == 0) return error.NoStructureFiles;
    try validateManifestPath(allocator, config, files);
    try std.Io.Dir.cwd().createDirPath(io, config.output_dir);
    try prepareRunMetadata(allocator, config);
    const results = try allocator.alloc(FileResult, files.len);
    defer allocator.free(results);
    @memset(results, .{});

    const cpu_count: u32 = @intCast(std.Thread.getCpuCount() catch 1);
    // A conservative default avoids multiplying whole-file decompression and
    // contact materialization by every logical CPU. Users can opt in to more.
    const requested = if (config.threads == 0) @min(cpu_count, 4) else @min(config.threads, cpu_count);
    const thread_count: u32 = @max(1, @min(requested, @as(u32, @intCast(files.len))));
    var context = WorkerContext{ .config = config, .files = files, .results = results };
    const threads = try allocator.alloc(std.Thread, thread_count);
    defer allocator.free(threads);
    var spawned: usize = 0;
    errdefer {
        _ = context.next.fetchAdd(files.len, .monotonic);
        for (threads[0..spawned]) |thread| thread.join();
    }
    for (threads) |*thread| {
        thread.* = try std.Thread.spawn(.{}, worker, .{&context});
        spawned += 1;
    }
    for (threads) |thread| thread.join();
    if (!config.quiet) std.debug.print("\n", .{});
    try writeManifest(allocator, config, files, results);

    var summary = Summary{ .total = files.len, .successful = 0, .skipped = 0, .failed = 0, .contacts = 0, .wall_time_ns = @intCast(started.untilNow(io, .awake).nanoseconds) };
    for (results) |result| switch (result.status) {
        .ok => {
            summary.successful += 1;
            summary.contacts += result.contact_count;
        },
        .skipped => {
            summary.skipped += 1;
            summary.contacts += result.contact_count;
        },
        .err => summary.failed += 1,
    };
    return summary;
}

test "structure filename filter" {
    try std.testing.expect(isStructureFile("x.pdb"));
    try std.testing.expect(isStructureFile("x.mmcif.gz"));
    try std.testing.expect(!isStructureFile("notes.tsv"));
}

test "parallel batch writes deterministic outputs and manifest" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.createDir(io, "in", @enumFromInt(0o755));
    const pdb =
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C  \n" ++
        "ATOM      2  CA  GLY A   2       3.000   0.000   0.000  1.00 10.00           C  \n";
    var input_dir = try tmp.dir.openDir(io, "in", .{});
    defer input_dir.close(io);
    try input_dir.writeFile(io, .{ .sub_path = "b.pdb", .data = pdb });
    try input_dir.writeFile(io, .{ .sub_path = "a.pdb", .data = pdb });
    const root_path = try tmp.dir.realPathFileAlloc(io, ".", allocator);
    defer allocator.free(root_path);
    const in_path = try std.fs.path.join(allocator, &.{ root_path, "in" });
    defer allocator.free(in_path);
    const out_path = try std.fs.path.join(allocator, &.{ root_path, "out" });
    defer allocator.free(out_path);
    const input_collision = try std.fs.path.join(allocator, &.{ in_path, "a.pdb" });
    defer allocator.free(input_collision);
    const collision_config = Config{ .input_dir = in_path, .output_dir = out_path, .manifest_path = input_collision, .quiet = true };
    try std.testing.expectError(error.ManifestPathCollision, run(allocator, &collision_config));
    const config = Config{ .input_dir = in_path, .output_dir = out_path, .threads = 2, .quiet = true };
    const summary = try run(allocator, &config);
    try std.testing.expectEqual(@as(usize, 2), summary.successful);
    try std.testing.expectEqual(@as(usize, 2), summary.contacts);
    const manifest_path = try std.fs.path.join(allocator, &.{ out_path, "manifest.tsv" });
    defer allocator.free(manifest_path);
    const manifest = try std.Io.Dir.cwd().readFileAlloc(io, manifest_path, allocator, .limited(4096));
    defer allocator.free(manifest);
    try std.testing.expect(std.mem.indexOf(u8, manifest, "a.pdb\ta.pdb.contacts.tsv\tok") != null);
    try std.testing.expect(std.mem.indexOf(u8, manifest, "b.pdb\tb.pdb.contacts.tsv\tok") != null);
    try std.testing.expect(std.mem.indexOf(u8, manifest, "a.pdb").? < std.mem.indexOf(u8, manifest, "b.pdb").?);

    const resume_config = Config{ .input_dir = in_path, .output_dir = out_path, .threads = 2, .quiet = true, .skip_existing = true };
    const resumed = try run(allocator, &resume_config);
    try std.testing.expectEqual(@as(usize, 2), resumed.skipped);
    try std.testing.expectEqual(@as(usize, 2), resumed.contacts);
    const resumed_manifest = try std.Io.Dir.cwd().readFileAlloc(io, manifest_path, allocator, .limited(4096));
    defer allocator.free(resumed_manifest);
    try std.testing.expect(std.mem.indexOf(u8, resumed_manifest, "a.pdb\ta.pdb.contacts.tsv\tskipped\t2\t2\t1") != null);

    // Existence alone is not sufficient for resume: a damaged output must be
    // recomputed, while the independently valid output remains skipped.
    const damaged_path = try std.fs.path.join(allocator, &.{ out_path, "a.pdb.contacts.tsv" });
    defer allocator.free(damaged_path);
    {
        const damaged = try std.Io.Dir.cwd().createFile(io, damaged_path, .{ .truncate = true });
        damaged.close(io);
    }
    const repaired = try run(allocator, &resume_config);
    try std.testing.expectEqual(@as(usize, 1), repaired.successful);
    try std.testing.expectEqual(@as(usize, 1), repaired.skipped);
    try std.testing.expectEqual(@as(usize, 2), repaired.contacts);

    // The marker also binds the output to the input content.
    try input_dir.writeFile(io, .{ .sub_path = "a.pdb", .data = pdb ++ "REMARK changed\n" });
    const changed = try run(allocator, &resume_config);
    try std.testing.expectEqual(@as(usize, 1), changed.successful);
    try std.testing.expectEqual(@as(usize, 1), changed.skipped);

    // A failed overwrite may replace run.tsv but must not make outputs from
    // the old scientific configuration eligible for a later resume.
    const missing_model_overwrite = Config{ .input_dir = in_path, .output_dir = out_path, .model = 2, .overwrite = true, .quiet = true };
    const overwrite_failed = try run(allocator, &missing_model_overwrite);
    try std.testing.expectEqual(@as(usize, 2), overwrite_failed.failed);
    const missing_model_resume = Config{ .input_dir = in_path, .output_dir = out_path, .model = 2, .skip_existing = true, .quiet = true };
    const resume_failed = try run(allocator, &missing_model_resume);
    try std.testing.expectEqual(@as(usize, 0), resume_failed.skipped);
    try std.testing.expectEqual(@as(usize, 2), resume_failed.failed);

    const incompatible = Config{ .input_dir = in_path, .output_dir = out_path, .cutoff = 5.0, .quiet = true, .skip_existing = true };
    try std.testing.expectError(error.ResumeConfigurationMismatch, run(allocator, &incompatible));
}
