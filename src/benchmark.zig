const std = @import("std");
const zcontact = @import("zcontact");

fn makeStructure(allocator: std.mem.Allocator, atom_count: usize) !zcontact.model.Structure {
    var structure = zcontact.model.Structure{};
    errdefer structure.deinit(allocator);
    try structure.atoms.ensureTotalCapacity(allocator, atom_count);
    for (0..atom_count) |i| {
        const fi: f64 = @floatFromInt(i);
        structure.atoms.appendAssumeCapacity(.{
            .serial = @intCast(i + 1),
            .name = try zcontact.model.Field.init("CA"),
            .residue_name = try zcontact.model.Field.init("GLY"),
            .chain = try zcontact.model.Field.init("A"),
            .residue_seq = @intCast(i + 1),
            .element = try zcontact.model.Field.init("C"),
            // A gently twisting chain gives bounded local density without a
            // best-case axis-aligned arrangement.
            .x = 2.0 * @cos(fi * 0.37),
            .y = 2.0 * @sin(fi * 0.37),
            .z = fi * 1.5,
            .residue_index = @intCast(i),
        });
    }
    return structure;
}

fn runCase(io: std.Io, allocator: std.mem.Allocator, atom_count: usize, iterations: usize) !void {
    var structure = try makeStructure(allocator, atom_count);
    defer structure.deinit(allocator);

    // One unmeasured warmup ensures lazy setup and code pages do not dominate.
    var warmup = try zcontact.contact.calculate(allocator, &structure, .atom, 4.0, "all", "all", .inter_residue);
    const contact_count = warmup.atom.items.len;
    warmup.deinit(allocator);

    const start = std.Io.Timestamp.now(io, .awake);
    for (0..iterations) |_| {
        var result = try zcontact.contact.calculate(allocator, &structure, .atom, 4.0, "all", "all", .inter_residue);
        std.mem.doNotOptimizeAway(result.atom.items.len);
        result.deinit(allocator);
    }
    const elapsed_ns: u64 = @intCast(start.untilNow(io, .awake).nanoseconds);
    const ms = @as(f64, @floatFromInt(elapsed_ns)) / @as(f64, @floatFromInt(iterations)) / 1_000_000.0;
    std.debug.print("{d}\t{d}\t{d:.3}\n", .{ atom_count, contact_count, ms });
}

pub fn main(init: std.process.Init) !void {
    const io = std.Io.Threaded.global_single_threaded.io();
    std.debug.print("atoms\tcontacts\tms_per_iteration\n", .{});
    try runCase(io, init.gpa, 1_000, 10);
    try runCase(io, init.gpa, 10_000, 5);
    try runCase(io, init.gpa, 100_000, 3);
}
