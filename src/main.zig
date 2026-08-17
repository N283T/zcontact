const std = @import("std");
const zcontact = @import("zcontact");
const build_options = @import("build_options");

const Config = struct {
    input: []const u8,
    output: ?[]const u8 = null,
    mode: zcontact.contact.Mode = .residue,
    cutoff: f64 = 4.0,
    select1: []const u8 = "polymer,heavy",
    select2: ?[]const u8 = null,
    scope: zcontact.contact.Scope = .inter_residue,
    model: u32 = 1,
    atoms_output: ?[]const u8 = null,
};

fn printStdout(comptime format: []const u8, args: anytype) void {
    const io = std.Io.Threaded.global_single_threaded.io();
    var buffer: [16 * 1024]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);
    writer.interface.print(format, args) catch std.process.exit(1);
    writer.interface.flush() catch std.process.exit(1);
}

fn usage(program: []const u8) void {
    printStdout(
        \\zcontact {s} — atomic and residue contacts in PDB/mmCIF structures
        \\
        \\USAGE:
        \\  {s} [OPTIONS] <structure.pdb|structure.cif>
        \\  {s} batch [OPTIONS] <input-directory> --output-dir DIR
        \\
        \\OPTIONS:
        \\  --mode atom|residue       Output atom pairs or residue pairs (default: residue)
        \\  --cutoff ANGSTROM         Inclusive Euclidean cutoff (default: 4.0)
        \\  --select EXPR             Symmetric selection (default: polymer,heavy)
        \\  --select1 EXPR            First side of an asymmetric selection
        \\  --select2 EXPR            Second side (default: same as select1)
        \\  --scope inter-residue|all Exclude/include contacts within one residue
        \\  --model N                 One 1-based PDB model number (default: 1)
        \\  -o, --output PATH         TSV destination (default: stdout)
        \\  --atoms-output PATH       Also write the selected atom inventory TSV
        \\  -h, --help                Show help
        \\  -V, --version             Show version
        \\
        \\SELECTION:
        \\  Comma means AND. Clauses: all, polymer, protein, hetero, ligand,
        \\  water, heavy, backbone, sidechain, ca, chain:ID, name:ATOM,
        \\  element:ELEMENT, resname:RESIDUE, resseq:N, icode:CODE.
        \\
        \\Examples:
        \\  {s} --mode residue --cutoff 4.0 structure.cif
        \\  {s} --mode atom --select1 chain:A --select2 chain:B structure.pdb
        \\  {s} batch structures/ -o contacts/ -j 8
        \\
    , .{ build_options.version, program, program, program, program, program });
}

fn batchUsage(program: []const u8) void {
    printStdout(
        \\USAGE:
        \\  {s} batch [OPTIONS] <input-directory> --output-dir DIR
        \\
        \\BATCH OPTIONS:
        \\  -o, --output-dir DIR       Per-structure TSV output directory (required)
        \\  --manifest PATH            Manifest TSV (default: DIR/manifest.tsv)
        \\  -j, --threads N            Worker threads (default: min(CPUs, 4))
        \\  --resume                   Reuse hash-validated matching outputs
        \\  --overwrite                Atomically replace existing outputs
        \\  --quiet                    Suppress progress and summary
        \\  --mode atom|residue        Contact output mode (default: residue)
        \\  --cutoff ANGSTROM          Inclusive cutoff (default: 4.0)
        \\  --select/--select1/--select2 EXPR
        \\  --scope inter-residue|all
        \\  --model N                  One model per structure (default: 1)
        \\
    , .{program});
}

fn needValue(args: []const []const u8, i: *usize, option: []const u8) []const u8 {
    i.* += 1;
    if (i.* >= args.len) {
        std.debug.print("Error: {s} requires a value\n", .{option});
        std.process.exit(2);
    }
    return args[i.*];
}

fn parseArgs(args: []const []const u8) ?Config {
    var input: ?[]const u8 = null;
    var config = Config{ .input = undefined };
    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            usage(args[0]);
            return null;
        }
        if (std.mem.eql(u8, arg, "-V") or std.mem.eql(u8, arg, "--version")) {
            printStdout("zcontact {s}\n", .{build_options.version});
            return null;
        }
        if (std.mem.eql(u8, arg, "--mode")) {
            const v = needValue(args, &i, arg);
            config.mode = if (std.mem.eql(u8, v, "atom")) .atom else if (std.mem.eql(u8, v, "residue")) .residue else {
                std.debug.print("Error: invalid --mode '{s}'\n", .{v});
                std.process.exit(2);
            };
        } else if (std.mem.eql(u8, arg, "--cutoff")) {
            const v = needValue(args, &i, arg);
            config.cutoff = std.fmt.parseFloat(f64, v) catch {
                std.debug.print("Error: invalid cutoff '{s}'\n", .{v});
                std.process.exit(2);
            };
            if (!std.math.isFinite(config.cutoff) or config.cutoff <= 0) {
                std.debug.print("Error: cutoff must be finite and > 0\n", .{});
                std.process.exit(2);
            }
        } else if (std.mem.eql(u8, arg, "--select")) {
            config.select1 = needValue(args, &i, arg);
            config.select2 = config.select1;
        } else if (std.mem.eql(u8, arg, "--select1")) {
            config.select1 = needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--select2")) {
            config.select2 = needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--scope")) {
            const v = needValue(args, &i, arg);
            config.scope = if (std.mem.eql(u8, v, "all")) .all else if (std.mem.eql(u8, v, "inter-residue")) .inter_residue else {
                std.debug.print("Error: invalid --scope '{s}'\n", .{v});
                std.process.exit(2);
            };
        } else if (std.mem.eql(u8, arg, "--model")) {
            const v = needValue(args, &i, arg);
            config.model = std.fmt.parseInt(u32, v, 10) catch {
                std.debug.print("Error: invalid model '{s}'\n", .{v});
                std.process.exit(2);
            };
            if (config.model == 0) {
                std.debug.print("Error: model must be >= 1\n", .{});
                std.process.exit(2);
            }
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            config.output = needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--atoms-output")) {
            config.atoms_output = needValue(args, &i, arg);
        } else if (std.mem.startsWith(u8, arg, "-")) {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(2);
        } else if (input == null) input = arg else {
            std.debug.print("Error: only one input structure is accepted\n", .{});
            std.process.exit(2);
        }
    }
    config.input = input orelse {
        usage(args[0]);
        std.process.exit(2);
    };
    zcontact.selection.validate(config.select1) catch {
        std.debug.print("Error: invalid selection '{s}'\n", .{config.select1});
        std.process.exit(2);
    };
    const s2 = config.select2 orelse config.select1;
    zcontact.selection.validate(s2) catch {
        std.debug.print("Error: invalid selection '{s}'\n", .{s2});
        std.process.exit(2);
    };
    return config;
}

fn parseBatchArgs(args: []const []const u8) ?zcontact.batch.Config {
    var input: ?[]const u8 = null;
    var output_dir: ?[]const u8 = null;
    var select2_explicit = false;
    var config = zcontact.batch.Config{ .input_dir = undefined, .output_dir = undefined };
    config.tool_version = build_options.version;
    var i: usize = 2;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            batchUsage(args[0]);
            return null;
        }
        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output-dir")) {
            output_dir = needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--manifest")) {
            config.manifest_path = needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--threads")) {
            const value = needValue(args, &i, arg);
            config.threads = std.fmt.parseInt(u32, value, 10) catch {
                std.debug.print("Error: invalid thread count '{s}'\n", .{value});
                std.process.exit(2);
            };
            if (config.threads == 0) {
                std.debug.print("Error: thread count must be >= 1\n", .{});
                std.process.exit(2);
            }
        } else if (std.mem.eql(u8, arg, "--resume")) {
            config.skip_existing = true;
        } else if (std.mem.eql(u8, arg, "--overwrite")) {
            config.overwrite = true;
        } else if (std.mem.eql(u8, arg, "--quiet")) {
            config.quiet = true;
        } else if (std.mem.eql(u8, arg, "--mode")) {
            const value = needValue(args, &i, arg);
            config.mode = if (std.mem.eql(u8, value, "atom")) .atom else if (std.mem.eql(u8, value, "residue")) .residue else {
                std.debug.print("Error: invalid --mode '{s}'\n", .{value});
                std.process.exit(2);
            };
        } else if (std.mem.eql(u8, arg, "--cutoff")) {
            const value = needValue(args, &i, arg);
            config.cutoff = std.fmt.parseFloat(f64, value) catch {
                std.debug.print("Error: invalid cutoff '{s}'\n", .{value});
                std.process.exit(2);
            };
            if (!std.math.isFinite(config.cutoff) or config.cutoff <= 0) {
                std.debug.print("Error: cutoff must be finite and > 0\n", .{});
                std.process.exit(2);
            }
        } else if (std.mem.eql(u8, arg, "--select")) {
            config.select1 = needValue(args, &i, arg);
            config.select2 = config.select1;
            select2_explicit = true;
        } else if (std.mem.eql(u8, arg, "--select1")) {
            config.select1 = needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--select2")) {
            config.select2 = needValue(args, &i, arg);
            select2_explicit = true;
        } else if (std.mem.eql(u8, arg, "--scope")) {
            const value = needValue(args, &i, arg);
            config.scope = if (std.mem.eql(u8, value, "all")) .all else if (std.mem.eql(u8, value, "inter-residue")) .inter_residue else {
                std.debug.print("Error: invalid --scope '{s}'\n", .{value});
                std.process.exit(2);
            };
        } else if (std.mem.eql(u8, arg, "--model")) {
            const value = needValue(args, &i, arg);
            config.model = std.fmt.parseInt(u32, value, 10) catch {
                std.debug.print("Error: invalid model '{s}'\n", .{value});
                std.process.exit(2);
            };
            if (config.model == 0) {
                std.debug.print("Error: model must be >= 1\n", .{});
                std.process.exit(2);
            }
        } else if (std.mem.startsWith(u8, arg, "-")) {
            std.debug.print("Error: unknown batch option '{s}'\n", .{arg});
            std.process.exit(2);
        } else if (input == null) input = arg else {
            std.debug.print("Error: batch accepts one input directory\n", .{});
            std.process.exit(2);
        }
    }
    config.input_dir = input orelse {
        batchUsage(args[0]);
        std.process.exit(2);
    };
    config.output_dir = output_dir orelse {
        std.debug.print("Error: batch requires --output-dir\n", .{});
        std.process.exit(2);
    };
    if (!select2_explicit) config.select2 = config.select1;
    if (config.skip_existing and config.overwrite) {
        std.debug.print("Error: --resume and --overwrite are mutually exclusive\n", .{});
        std.process.exit(2);
    }
    zcontact.selection.validate(config.select1) catch {
        std.debug.print("Error: invalid selection '{s}'\n", .{config.select1});
        std.process.exit(2);
    };
    zcontact.selection.validate(config.select2) catch {
        std.debug.print("Error: invalid selection '{s}'\n", .{config.select2});
        std.process.exit(2);
    };
    return config;
}

fn temporaryPath(allocator: std.mem.Allocator, io: std.Io, destination: []const u8) ![]u8 {
    const absolute = try std.fs.path.resolve(allocator, &.{destination});
    defer allocator.free(absolute);
    return std.fmt.allocPrint(allocator, "{s}.tmp-{d}-{d}", .{ absolute, std.Thread.getCurrentId(), std.Io.Timestamp.now(io, .awake).nanoseconds });
}

fn canonicalPath(allocator: std.mem.Allocator, io: std.Io, path: []const u8) ![]u8 {
    const real = std.Io.Dir.cwd().realPathFileAlloc(io, path, allocator) catch |err| switch (err) {
        error.FileNotFound => return std.fs.path.resolve(allocator, &.{path}),
        else => return err,
    };
    defer allocator.free(real);
    return allocator.dupe(u8, real);
}

fn writeAtomsAtomic(allocator: std.mem.Allocator, path: []const u8, structure: *const zcontact.model.Structure, select1: []const u8, select2: []const u8) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const temp = try temporaryPath(allocator, io, path);
    defer allocator.free(temp);
    errdefer std.Io.Dir.cwd().deleteFile(io, temp) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, temp, .{ .exclusive = true });
        defer file.close(io);
        var buffer: [64 * 1024]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try zcontact.output.writeSelectedAtoms(&writer.interface, structure, select1, select2);
        try writer.interface.flush();
    }
    const destination = try std.fs.path.resolve(allocator, &.{path});
    defer allocator.free(destination);
    try std.Io.Dir.renameAbsolute(temp, destination, io);
}

fn writeContactsAtomic(allocator: std.mem.Allocator, path: []const u8, structure: *const zcontact.model.Structure, result: *const zcontact.contact.Result) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const temp = try temporaryPath(allocator, io, path);
    defer allocator.free(temp);
    errdefer std.Io.Dir.cwd().deleteFile(io, temp) catch {};
    {
        const file = try std.Io.Dir.cwd().createFile(io, temp, .{ .exclusive = true });
        defer file.close(io);
        var buffer: [64 * 1024]u8 = undefined;
        var writer = file.writer(io, &buffer);
        try zcontact.output.writeResult(&writer.interface, structure, result);
        try writer.interface.flush();
    }
    const destination = try std.fs.path.resolve(allocator, &.{path});
    defer allocator.free(destination);
    try std.Io.Dir.renameAbsolute(temp, destination, io);
}

fn run(allocator: std.mem.Allocator, config: Config) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const input_resolved = try canonicalPath(allocator, io, config.input);
    defer allocator.free(input_resolved);
    if (config.output) |path| {
        const resolved = try canonicalPath(allocator, io, path);
        defer allocator.free(resolved);
        if (std.mem.eql(u8, input_resolved, resolved)) return error.PathCollision;
    }
    if (config.atoms_output) |path| {
        const resolved = try canonicalPath(allocator, io, path);
        defer allocator.free(resolved);
        if (std.mem.eql(u8, input_resolved, resolved)) return error.PathCollision;
        if (config.output) |output_path| {
            const output_resolved = try canonicalPath(allocator, io, output_path);
            defer allocator.free(output_resolved);
            if (std.mem.eql(u8, output_resolved, resolved)) return error.PathCollision;
        }
    }
    const source = try zcontact.gzip.readFile(allocator, config.input, zcontact.gzip.default_max_size);
    defer allocator.free(source);
    const format = try zcontact.parser.detectFormat(config.input, source);
    var structure = try zcontact.parser.parse(allocator, source, format, config.model);
    defer structure.deinit(allocator);
    var result = try zcontact.contact.calculate(allocator, &structure, config.mode, config.cutoff, config.select1, config.select2 orelse config.select1, config.scope);
    defer result.deinit(allocator);

    if (config.atoms_output) |path|
        try writeAtomsAtomic(allocator, path, &structure, config.select1, config.select2 orelse config.select1);

    var buffer: [64 * 1024]u8 = undefined;
    if (config.output) |path| {
        try writeContactsAtomic(allocator, path, &structure, &result);
    } else {
        const stdout = std.Io.File.stdout();
        var sw = stdout.writer(io, &buffer);
        try zcontact.output.writeResult(&sw.interface, &structure, &result);
        try sw.interface.flush();
    }
}

pub fn main(init: std.process.Init) void {
    const args = init.minimal.args.toSlice(init.arena.allocator()) catch {
        std.debug.print("Error: cannot read arguments\n", .{});
        std.process.exit(2);
    };
    if (args.len > 1 and std.mem.eql(u8, args[1], "batch")) {
        const batch_config = parseBatchArgs(args) orelse return;
        const summary = zcontact.batch.run(init.gpa, &batch_config) catch |err| {
            std.debug.print("Error: batch: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
        if (!batch_config.quiet) std.debug.print("Batch: {d} total, {d} ok, {d} skipped, {d} failed, {d} contacts, {d:.3} s\n", .{ summary.total, summary.successful, summary.skipped, summary.failed, summary.contacts, @as(f64, @floatFromInt(summary.wall_time_ns)) / 1_000_000_000.0 });
        if (summary.failed != 0) std.process.exit(1);
        return;
    }
    const config = parseArgs(args) orelse return;
    run(init.gpa, config) catch |err| {
        std.debug.print("Error: {s}: {s}\n", .{ config.input, @errorName(err) });
        std.process.exit(1);
    };
}
