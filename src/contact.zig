const std = @import("std");
const model = @import("model.zig");
const selection = @import("selection.zig");

pub const Scope = enum { inter_residue, all };
pub const Mode = enum { atom, residue };

pub const AtomContact = struct { a: u32, b: u32, distance: f64 };
pub const ResidueContact = struct {
    atom_a: u32,
    atom_b: u32,
    residue_a: u32,
    residue_b: u32,
    distance: f64,
};

pub const Result = union(Mode) {
    atom: std.ArrayListUnmanaged(AtomContact),
    residue: std.ArrayListUnmanaged(ResidueContact),

    pub fn deinit(self: *Result, allocator: std.mem.Allocator) void {
        switch (self.*) {
            .atom => |*v| v.deinit(allocator),
            .residue => |*v| v.deinit(allocator),
        }
    }
};

const ResiduePair = struct { a: u32, b: u32 };
const Cell = struct { x: i64, y: i64, z: i64 };
const CellEntry = struct { cell: Cell, atom: u32 };
const CellRange = struct { start: usize, end: usize };

pub fn calculate(allocator: std.mem.Allocator, structure: *const model.Structure, mode: Mode, cutoff: f64, select1: []const u8, select2: []const u8, scope: Scope) !Result {
    if (!std.math.isFinite(cutoff) or cutoff <= 0 or !std.math.isFinite(cutoff * cutoff)) return error.InvalidCutoff;
    try selection.validate(select1);
    try selection.validate(select2);
    const atoms = structure.atoms.items;
    const cutoff2 = cutoff * cutoff;
    const mask1 = try allocator.alloc(bool, atoms.len);
    defer allocator.free(mask1);
    const mask2 = try allocator.alloc(bool, atoms.len);
    defer allocator.free(mask2);
    for (atoms, 0..) |atom, i| {
        mask1[i] = selection.matches(atom, select1);
        mask2[i] = selection.matches(atom, select2);
    }

    // Sort atoms into cutoff-sized cells, then visit the 27 neighboring cells.
    // This preserves exact cutoff semantics while avoiding a quadratic scan.
    const entries = try allocator.alloc(CellEntry, atoms.len);
    defer allocator.free(entries);
    for (atoms, 0..) |atom, i| entries[i] = .{ .cell = try cellFor(atom, cutoff), .atom = @intCast(i) };
    std.mem.sort(CellEntry, entries, {}, cellEntryLessThan);

    var ranges = std.AutoHashMapUnmanaged(Cell, CellRange).empty;
    defer ranges.deinit(allocator);
    var start: usize = 0;
    while (start < entries.len) {
        var end = start + 1;
        while (end < entries.len and cellEqual(entries[start].cell, entries[end].cell)) end += 1;
        try ranges.put(allocator, entries[start].cell, .{ .start = start, .end = end });
        start = end;
    }

    if (mode == .atom) {
        var contacts = std.ArrayListUnmanaged(AtomContact).empty;
        errdefer contacts.deinit(allocator);
        for (atoms, 0..) |a, i| {
            const home = try cellFor(a, cutoff);
            for ([_]i64{ -1, 0, 1 }) |dx| for ([_]i64{ -1, 0, 1 }) |dy| for ([_]i64{ -1, 0, 1 }) |dz| {
                const range = ranges.get(.{ .x = home.x + dx, .y = home.y + dy, .z = home.z + dz }) orelse continue;
                for (entries[range.start..range.end]) |entry| {
                    const j: usize = entry.atom;
                    if (j <= i) continue;
                    const b = atoms[j];
                    const oriented = orientPair(@intCast(i), @intCast(j), mask1[i], mask2[i], mask1[j], mask2[j]) orelse continue;
                    if (scope == .inter_residue and a.residue_index == b.residue_index) continue;
                    const d2 = distanceSquared(a, b);
                    if (d2 <= cutoff2) try contacts.append(allocator, .{ .a = oriented.a, .b = oriented.b, .distance = @sqrt(d2) });
                }
            };
        }
        std.mem.sort(AtomContact, contacts.items, {}, atomContactLessThan);
        return .{ .atom = contacts };
    }

    var index = std.AutoHashMapUnmanaged(ResiduePair, usize).empty;
    defer index.deinit(allocator);
    var contacts = std.ArrayListUnmanaged(ResidueContact).empty;
    errdefer contacts.deinit(allocator);
    for (atoms, 0..) |a, i| {
        const home = try cellFor(a, cutoff);
        for ([_]i64{ -1, 0, 1 }) |dx| for ([_]i64{ -1, 0, 1 }) |dy| for ([_]i64{ -1, 0, 1 }) |dz| {
            const range = ranges.get(.{ .x = home.x + dx, .y = home.y + dy, .z = home.z + dz }) orelse continue;
            for (entries[range.start..range.end]) |entry| {
                const j: usize = entry.atom;
                if (j <= i) continue;
                const b = atoms[j];
                const oriented = orientPair(@intCast(i), @intCast(j), mask1[i], mask2[i], mask1[j], mask2[j]) orelse continue;
                if (scope == .inter_residue and a.residue_index == b.residue_index) continue;
                const d2 = distanceSquared(a, b);
                if (d2 > cutoff2) continue;
                const key = ResiduePair{ .a = @min(a.residue_index, b.residue_index), .b = @max(a.residue_index, b.residue_index) };
                const d = @sqrt(d2);
                if (index.get(key)) |idx| {
                    const old = contacts.items[idx];
                    if (d < old.distance or (d == old.distance and atomPairLessThan(oriented.a, oriented.b, old.atom_a, old.atom_b)))
                        contacts.items[idx] = .{ .atom_a = oriented.a, .atom_b = oriented.b, .residue_a = key.a, .residue_b = key.b, .distance = d };
                } else {
                    try index.put(allocator, key, contacts.items.len);
                    try contacts.append(allocator, .{ .atom_a = oriented.a, .atom_b = oriented.b, .residue_a = key.a, .residue_b = key.b, .distance = d });
                }
            }
        };
    }
    std.mem.sort(ResidueContact, contacts.items, {}, residueContactLessThan);
    return .{ .residue = contacts };
}

fn cellFor(atom: model.Atom, cutoff: f64) !Cell {
    return .{ .x = try cellCoordinate(atom.x, cutoff), .y = try cellCoordinate(atom.y, cutoff), .z = try cellCoordinate(atom.z, cutoff) };
}

fn cellCoordinate(value: f64, cutoff: f64) !i64 {
    const scaled = @floor(value / cutoff);
    // Reserve one cell at each extreme so neighbor offsets cannot overflow.
    const lower: f64 = @floatFromInt(std.math.minInt(i64) + 1);
    const upper: f64 = @floatFromInt(std.math.maxInt(i64) - 1);
    if (!std.math.isFinite(scaled) or scaled < lower or scaled > upper) return error.CoordinateOutOfGrid;
    return @intFromFloat(scaled);
}

fn cellEqual(a: Cell, b: Cell) bool {
    return a.x == b.x and a.y == b.y and a.z == b.z;
}

fn cellEntryLessThan(_: void, a: CellEntry, b: CellEntry) bool {
    if (a.cell.x != b.cell.x) return a.cell.x < b.cell.x;
    if (a.cell.y != b.cell.y) return a.cell.y < b.cell.y;
    if (a.cell.z != b.cell.z) return a.cell.z < b.cell.z;
    return a.atom < b.atom;
}

fn atomPairLessThan(a1: u32, b1: u32, a2: u32, b2: u32) bool {
    return a1 < a2 or (a1 == a2 and b1 < b2);
}

fn atomContactLessThan(_: void, a: AtomContact, b: AtomContact) bool {
    return atomPairLessThan(a.a, a.b, b.a, b.b);
}

fn residueContactLessThan(_: void, a: ResidueContact, b: ResidueContact) bool {
    return a.residue_a < b.residue_a or (a.residue_a == b.residue_a and a.residue_b < b.residue_b);
}

const OrientedPair = struct { a: u32, b: u32 };

fn orientPair(a: u32, b: u32, a1: bool, a2: bool, b1: bool, b2: bool) ?OrientedPair {
    if (a1 and b2) return .{ .a = a, .b = b };
    if (b1 and a2) return .{ .a = b, .b = a };
    return null;
}

fn distanceSquared(a: model.Atom, b: model.Atom) f64 {
    const dx = a.x - b.x;
    const dy = a.y - b.y;
    const dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

test "atom and residue minimum-distance contacts" {
    var s = model.Structure{};
    defer s.deinit(std.testing.allocator);
    for ([_]f64{ 0, 1, 3 }) |x| try s.atoms.append(std.testing.allocator, .{ .serial = @intFromFloat(x + 1), .name = try model.Field.init("CA"), .residue_name = try model.Field.init("ALA"), .chain = try model.Field.init("A"), .residue_seq = if (x < 2) 1 else 2, .element = try model.Field.init("C"), .x = x, .y = 0, .z = 0, .residue_index = if (x < 2) 0 else 1 });
    var atoms = try calculate(std.testing.allocator, &s, .atom, 2.1, "all", "all", .inter_residue);
    defer atoms.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), atoms.atom.items.len);
    var residues = try calculate(std.testing.allocator, &s, .residue, 4, "all", "all", .inter_residue);
    defer residues.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), residues.residue.items.len);
    try std.testing.expectApproxEqAbs(@as(f64, 2), residues.residue.items[0].distance, 1e-9);
}

test "cell list matches brute force across negative and boundary cells" {
    const coordinates = [_][3]f64{
        .{ -4.0, 0.0, 0.0 },   .{ -0.01, 0.0, 0.0 },  .{ 0.0, 0.0, 0.0 },
        .{ 3.99, 0.0, 0.0 },   .{ 4.0, 0.0, 0.0 },    .{ 7.0, 2.0, 1.0 },
        .{ 20.0, 20.0, 20.0 }, .{ -3.0, -2.0, -1.0 },
    };
    var s = model.Structure{};
    defer s.deinit(std.testing.allocator);
    for (coordinates, 0..) |xyz, i| try s.atoms.append(std.testing.allocator, .{
        .serial = @intCast(i + 1),
        .name = try model.Field.init("CA"),
        .residue_name = try model.Field.init("GLY"),
        .chain = try model.Field.init("A"),
        .residue_seq = @intCast(i + 1),
        .element = try model.Field.init("C"),
        .x = xyz[0],
        .y = xyz[1],
        .z = xyz[2],
        .residue_index = @intCast(i),
    });

    var result = try calculate(std.testing.allocator, &s, .atom, 4.0, "all", "all", .inter_residue);
    defer result.deinit(std.testing.allocator);
    var expected_count: usize = 0;
    var actual_index: usize = 0;
    for (s.atoms.items, 0..) |a, i| for (s.atoms.items[i + 1 ..], i + 1..) |b, j| {
        const d2 = distanceSquared(a, b);
        if (d2 <= 16.0) {
            expected_count += 1;
            const actual = result.atom.items[actual_index];
            try std.testing.expectEqual(@as(u32, @intCast(i)), actual.a);
            try std.testing.expectEqual(@as(u32, @intCast(j)), actual.b);
            try std.testing.expectApproxEqAbs(@sqrt(d2), actual.distance, 1e-12);
            actual_index += 1;
        }
    };
    try std.testing.expectEqual(expected_count, result.atom.items.len);
}

test "library rejects invalid cutoff and selection" {
    const empty = model.Structure{};
    try std.testing.expectError(error.InvalidCutoff, calculate(std.testing.allocator, &empty, .atom, 0, "all", "all", .all));
    try std.testing.expectError(error.InvalidSelection, calculate(std.testing.allocator, &empty, .atom, 4, "unknown", "all", .all));
}

test "asymmetric selection orients output as select1 then select2" {
    var s = model.Structure{};
    defer s.deinit(std.testing.allocator);
    try s.atoms.append(std.testing.allocator, .{ .serial = 1, .name = try model.Field.init("CA"), .residue_name = try model.Field.init("GLY"), .chain = try model.Field.init("B"), .residue_seq = 1, .element = try model.Field.init("C"), .x = 0, .y = 0, .z = 0, .residue_index = 0 });
    try s.atoms.append(std.testing.allocator, .{ .serial = 2, .name = try model.Field.init("CA"), .residue_name = try model.Field.init("GLY"), .chain = try model.Field.init("A"), .residue_seq = 2, .element = try model.Field.init("C"), .x = 3, .y = 0, .z = 0, .residue_index = 1 });
    var result = try calculate(std.testing.allocator, &s, .atom, 4, "chain:A", "chain:B", .inter_residue);
    defer result.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(u32, 1), result.atom.items[0].a);
    try std.testing.expectEqual(@as(u32, 0), result.atom.items[0].b);
}
