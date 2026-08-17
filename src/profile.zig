//! Development-only, real-file stage profiler. This is intentionally separate
//! from the stable `zcontact` CLI; invoke it with `zig build profile -- ...`.

const std = @import("std");
const zcontact = @import("zcontact");

const Config = struct {
    mode: zcontact.contact.Mode = .residue,
    cutoff: f64 = 4.0,
    select1: []const u8 = "polymer,heavy",
    select2: ?[]const u8 = null,
    scope: zcontact.contact.Scope = .inter_residue,
    model: u32 = 1,
    iterations: usize = 5,
};

fn elapsedNs(start: std.Io.Timestamp, io: std.Io) u64 {
    return @intCast(start.untilNow(io, .awake).nanoseconds);
}

fn needValue(args: []const []const u8, i: *usize, option: []const u8) ![]const u8 {
    i.* += 1;
    if (i.* >= args.len) {
        std.debug.print("Error: {s} requires a value\n", .{option});
        return error.InvalidArguments;
    }
    return args[i.*];
}

fn printUsage(program: []const u8) void {
    std.debug.print(
        \\USAGE: {s} [OPTIONS] <structure> [structure ...]
        \\
        \\Development-only ReleaseFast stage profiler. Each file gets one
        \\unmeasured warmup followed by TSV rows for measured iterations.
        \\
        \\OPTIONS:
        \\  --mode atom|residue       Contact mode (default: residue)
        \\  --cutoff ANGSTROM         Inclusive cutoff (default: 4.0)
        \\  --select EXPR             Symmetric selection
        \\  --select1/--select2 EXPR  Asymmetric selections
        \\  --scope inter-residue|all Contact scope
        \\  --model N                 Model number (default: 1)
        \\  --iterations N            Measured runs per file (default: 5)
        \\  -h, --help                Show this help
        \\
    , .{program});
}

fn parseOptions(args: []const []const u8, first_path: *usize) !?Config {
    var config = Config{};
    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (!std.mem.startsWith(u8, arg, "-")) break;
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printUsage(args[0]);
            return null;
        } else if (std.mem.eql(u8, arg, "--mode")) {
            const value = try needValue(args, &i, arg);
            config.mode = if (std.mem.eql(u8, value, "atom")) .atom else if (std.mem.eql(u8, value, "residue")) .residue else return error.InvalidArguments;
        } else if (std.mem.eql(u8, arg, "--cutoff")) {
            config.cutoff = try std.fmt.parseFloat(f64, try needValue(args, &i, arg));
            if (!std.math.isFinite(config.cutoff) or config.cutoff <= 0) return error.InvalidArguments;
        } else if (std.mem.eql(u8, arg, "--select")) {
            config.select1 = try needValue(args, &i, arg);
            config.select2 = config.select1;
        } else if (std.mem.eql(u8, arg, "--select1")) {
            config.select1 = try needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--select2")) {
            config.select2 = try needValue(args, &i, arg);
        } else if (std.mem.eql(u8, arg, "--scope")) {
            const value = try needValue(args, &i, arg);
            config.scope = if (std.mem.eql(u8, value, "all")) .all else if (std.mem.eql(u8, value, "inter-residue")) .inter_residue else return error.InvalidArguments;
        } else if (std.mem.eql(u8, arg, "--model")) {
            config.model = try std.fmt.parseInt(u32, try needValue(args, &i, arg), 10);
            if (config.model == 0) return error.InvalidArguments;
        } else if (std.mem.eql(u8, arg, "--iterations")) {
            config.iterations = try std.fmt.parseInt(usize, try needValue(args, &i, arg), 10);
            if (config.iterations == 0) return error.InvalidArguments;
        } else {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            return error.InvalidArguments;
        }
    }
    try zcontact.selection.validate(config.select1);
    try zcontact.selection.validate(config.select2 orelse config.select1);
    first_path.* = i;
    return config;
}

fn runOnce(allocator: std.mem.Allocator, path: []const u8, config: Config, iteration: usize, emit: bool, writer: *std.Io.Writer) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    const pipeline_start = std.Io.Timestamp.now(io, .awake);

    var start = std.Io.Timestamp.now(io, .awake);
    const source = try zcontact.gzip.readFile(allocator, path, zcontact.gzip.default_max_size);
    defer allocator.free(source);
    const read_ns = elapsedNs(start, io);

    start = std.Io.Timestamp.now(io, .awake);
    const format = try zcontact.parser.detectFormat(path, source);
    var parser_profile = zcontact.parser.Profile{};
    var structure = try zcontact.parser.parseProfiled(allocator, source, format, config.model, &parser_profile);
    defer structure.deinit(allocator);
    const parse_ns = elapsedNs(start, io);

    var contact_profile = zcontact.contact.Profile{};
    var result = try zcontact.contact.calculateProfiled(allocator, &structure, config.mode, config.cutoff, config.select1, config.select2 orelse config.select1, config.scope, &contact_profile);
    defer result.deinit(allocator);

    start = std.Io.Timestamp.now(io, .awake);
    var discard_buffer: [64 * 1024]u8 = undefined;
    var discarding: std.Io.Writer.Discarding = .init(&discard_buffer);
    try zcontact.output.writeResult(&discarding.writer, &structure, &result);
    try discarding.writer.flush();
    const output_ns = elapsedNs(start, io);
    const total_ns = elapsedNs(pipeline_start, io);

    if (!emit) return;
    const input_file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer input_file.close(io);
    const input_bytes = (try input_file.stat(io)).size;
    try writer.print("{s}\t{s}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\t{d}\n", .{
        path,
        @tagName(config.mode),
        iteration,
        input_bytes,
        source.len,
        discarding.fullCount(),
        total_ns,
        read_ns,
        parse_ns,
        parser_profile.raw_parse_ns,
        parser_profile.altloc_ns,
        parser_profile.residue_assignment_ns,
        contact_profile.selection_ns,
        contact_profile.grid_ns,
        contact_profile.search_ns,
        contact_profile.sort_ns,
        output_ns,
        contact_profile.atoms,
        contact_profile.selected1_atoms,
        contact_profile.selected2_atoms,
        contact_profile.occupied_cells,
        contact_profile.candidate_pairs,
        contact_profile.selection_pairs,
        contact_profile.distance_evaluations,
        contact_profile.accepted_atom_pairs,
        contact_profile.result_count,
    });
}

pub fn main(init: std.process.Init) !void {
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    var first_path: usize = 0;
    const config = (try parseOptions(args, &first_path)) orelse return;
    if (first_path >= args.len) {
        printUsage(args[0]);
        return error.InvalidArguments;
    }
    const io = std.Io.Threaded.global_single_threaded.io();
    var stdout_buffer: [64 * 1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(io, &stdout_buffer);
    try stdout.interface.writeAll("path\tmode\titeration\tinput_bytes\tsource_bytes\toutput_bytes\ttotal_ns\tread_decompress_ns\tparse_ns\traw_parse_ns\taltloc_ns\tresidue_assignment_ns\tselection_ns\tgrid_ns\tsearch_aggregate_ns\tsort_ns\toutput_ns\tatoms\tselected1_atoms\tselected2_atoms\toccupied_cells\tcandidate_pairs\tselection_pairs\tdistance_evaluations\taccepted_atom_pairs\tresults\n");
    for (args[first_path..]) |path| {
        try runOnce(init.gpa, path, config, 0, false, &stdout.interface);
        for (1..config.iterations + 1) |iteration| try runOnce(init.gpa, path, config, iteration, true, &stdout.interface);
    }
    try stdout.interface.flush();
}
