const std = @import("std");
const model = @import("model.zig");
const selection = @import("selection.zig");

pub const Scope = enum { inter_residue, all };
pub const Mode = enum { atom, residue };

/// Optional development-time instrumentation for the deterministic cell-list
/// engine. `calculate` does not collect these counters; profiling tools must
/// opt in through `calculateProfiled`.
pub const Profile = struct {
    selection_ns: u64 = 0,
    grid_ns: u64 = 0,
    search_ns: u64 = 0,
    sort_ns: u64 = 0,
    atoms: usize = 0,
    selected1_atoms: usize = 0,
    selected2_atoms: usize = 0,
    occupied_cells: usize = 0,
    candidate_pairs: u64 = 0,
    selection_pairs: u64 = 0,
    distance_evaluations: u64 = 0,
    accepted_atom_pairs: u64 = 0,
    result_count: usize = 0,
};

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
    return calculateImpl(allocator, structure, mode, cutoff, select1, select2, scope, null);
}

pub fn calculateProfiled(allocator: std.mem.Allocator, structure: *const model.Structure, mode: Mode, cutoff: f64, select1: []const u8, select2: []const u8, scope: Scope, profile: *Profile) !Result {
    profile.* = .{};
    return calculateImpl(allocator, structure, mode, cutoff, select1, select2, scope, profile);
}

fn calculateImpl(allocator: std.mem.Allocator, structure: *const model.Structure, mode: Mode, cutoff: f64, select1: []const u8, select2: []const u8, scope: Scope, profile: ?*Profile) !Result {
    if (!std.math.isFinite(cutoff) or cutoff <= 0 or !std.math.isFinite(cutoff * cutoff)) return error.InvalidCutoff;
    try selection.validate(select1);
    try selection.validate(select2);
    const atoms = structure.atoms.items;
    if (profile) |p| p.atoms = atoms.len;
    const io = std.Io.Threaded.global_single_threaded.io();
    var stage_start: std.Io.Timestamp = undefined;
    if (profile != null) stage_start = std.Io.Timestamp.now(io, .awake);
    const cutoff2 = cutoff * cutoff;
    const mask1 = try allocator.alloc(bool, atoms.len);
    defer allocator.free(mask1);
    const mask2 = try allocator.alloc(bool, atoms.len);
    defer allocator.free(mask2);
    for (atoms, 0..) |atom, i| {
        mask1[i] = selection.matches(atom, select1);
        mask2[i] = selection.matches(atom, select2);
        if (profile) |p| {
            if (mask1[i]) p.selected1_atoms += 1;
            if (mask2[i]) p.selected2_atoms += 1;
        }
    }
    if (profile) |p| {
        p.selection_ns = elapsedNs(stage_start, io);
        stage_start = std.Io.Timestamp.now(io, .awake);
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
    if (profile) |p| {
        p.occupied_cells = ranges.count();
        p.grid_ns = elapsedNs(stage_start, io);
        stage_start = std.Io.Timestamp.now(io, .awake);
    }

    if (mode == .atom) {
        var contacts = std.ArrayListUnmanaged(AtomContact).empty;
        errdefer contacts.deinit(allocator);
        start = 0;
        while (start < entries.len) {
            const home = entries[start].cell;
            const home_range = ranges.get(home).?;
            start = home_range.end;
            for ([_]i64{ -1, 0, 1 }) |dx| for ([_]i64{ -1, 0, 1 }) |dy| for ([_]i64{ -1, 0, 1 }) |dz| {
                const neighbor = Cell{ .x = home.x + dx, .y = home.y + dy, .z = home.z + dz };
                if (cellLessThan(neighbor, home)) continue;
                const range = ranges.get(neighbor) orelse continue;
                for (entries[home_range.start..home_range.end]) |home_entry| for (entries[range.start..range.end]) |entry| {
                    if (cellEqual(home, neighbor) and entry.atom <= home_entry.atom) continue;
                    const i: usize = @min(home_entry.atom, entry.atom);
                    const j: usize = @max(home_entry.atom, entry.atom);
                    if (profile) |p| p.candidate_pairs += 1;
                    const a = atoms[i];
                    const b = atoms[j];
                    const oriented = orientPair(@intCast(i), @intCast(j), mask1[i], mask2[i], mask1[j], mask2[j]) orelse continue;
                    if (profile) |p| p.selection_pairs += 1;
                    if (scope == .inter_residue and a.residue_index == b.residue_index) continue;
                    if (profile) |p| p.distance_evaluations += 1;
                    const d2 = distanceSquared(a, b);
                    if (d2 <= cutoff2) {
                        if (profile) |p| p.accepted_atom_pairs += 1;
                        try contacts.append(allocator, .{ .a = oriented.a, .b = oriented.b, .distance = @sqrt(d2) });
                    }
                };
            };
        }
        if (profile) |p| {
            p.search_ns = elapsedNs(stage_start, io);
            stage_start = std.Io.Timestamp.now(io, .awake);
        }
        std.mem.sort(AtomContact, contacts.items, {}, atomContactLessThan);
        if (profile) |p| {
            p.sort_ns = elapsedNs(stage_start, io);
            p.result_count = contacts.items.len;
        }
        return .{ .atom = contacts };
    }

    var index = std.AutoHashMapUnmanaged(ResiduePair, usize).empty;
    defer index.deinit(allocator);
    var contacts = std.ArrayListUnmanaged(ResidueContact).empty;
    errdefer contacts.deinit(allocator);
    start = 0;
    while (start < entries.len) {
        const home = entries[start].cell;
        const home_range = ranges.get(home).?;
        start = home_range.end;
        for ([_]i64{ -1, 0, 1 }) |dx| for ([_]i64{ -1, 0, 1 }) |dy| for ([_]i64{ -1, 0, 1 }) |dz| {
            const neighbor = Cell{ .x = home.x + dx, .y = home.y + dy, .z = home.z + dz };
            if (cellLessThan(neighbor, home)) continue;
            const range = ranges.get(neighbor) orelse continue;
            for (entries[home_range.start..home_range.end]) |home_entry| for (entries[range.start..range.end]) |entry| {
                if (cellEqual(home, neighbor) and entry.atom <= home_entry.atom) continue;
                const i: usize = @min(home_entry.atom, entry.atom);
                const j: usize = @max(home_entry.atom, entry.atom);
                if (profile) |p| p.candidate_pairs += 1;
                const a = atoms[i];
                const b = atoms[j];
                const oriented = orientPair(@intCast(i), @intCast(j), mask1[i], mask2[i], mask1[j], mask2[j]) orelse continue;
                if (profile) |p| p.selection_pairs += 1;
                if (scope == .inter_residue and a.residue_index == b.residue_index) continue;
                if (profile) |p| p.distance_evaluations += 1;
                const d2 = distanceSquared(a, b);
                if (d2 > cutoff2) continue;
                if (profile) |p| p.accepted_atom_pairs += 1;
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
            };
        };
    }
    if (profile) |p| {
        p.search_ns = elapsedNs(stage_start, io);
        stage_start = std.Io.Timestamp.now(io, .awake);
    }
    std.mem.sort(ResidueContact, contacts.items, {}, residueContactLessThan);
    if (profile) |p| {
        p.sort_ns = elapsedNs(stage_start, io);
        p.result_count = contacts.items.len;
    }
    return .{ .residue = contacts };
}

fn elapsedNs(start: std.Io.Timestamp, io: std.Io) u64 {
    return @intCast(start.untilNow(io, .awake).nanoseconds);
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

fn cellLessThan(a: Cell, b: Cell) bool {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.z < b.z;
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

fn calculateBruteForce(allocator: std.mem.Allocator, structure: *const model.Structure, mode: Mode, cutoff: f64, select1: []const u8, select2: []const u8, scope: Scope) !Result {
    const atoms = structure.atoms.items;
    const cutoff2 = cutoff * cutoff;
    if (mode == .atom) {
        var contacts = std.ArrayListUnmanaged(AtomContact).empty;
        errdefer contacts.deinit(allocator);
        for (atoms, 0..) |a, i| for (atoms[i + 1 ..], i + 1..) |b, j| {
            const oriented = orientPair(@intCast(i), @intCast(j), selection.matches(a, select1), selection.matches(a, select2), selection.matches(b, select1), selection.matches(b, select2)) orelse continue;
            if (scope == .inter_residue and a.residue_index == b.residue_index) continue;
            const d2 = distanceSquared(a, b);
            if (d2 <= cutoff2) try contacts.append(allocator, .{ .a = oriented.a, .b = oriented.b, .distance = @sqrt(d2) });
        };
        std.mem.sort(AtomContact, contacts.items, {}, atomContactLessThan);
        return .{ .atom = contacts };
    }

    var index = std.AutoHashMapUnmanaged(ResiduePair, usize).empty;
    defer index.deinit(allocator);
    var contacts = std.ArrayListUnmanaged(ResidueContact).empty;
    errdefer contacts.deinit(allocator);
    for (atoms, 0..) |a, i| for (atoms[i + 1 ..], i + 1..) |b, j| {
        const oriented = orientPair(@intCast(i), @intCast(j), selection.matches(a, select1), selection.matches(a, select2), selection.matches(b, select1), selection.matches(b, select2)) orelse continue;
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
    };
    std.mem.sort(ResidueContact, contacts.items, {}, residueContactLessThan);
    return .{ .residue = contacts };
}

fn expectSameResult(expected: *const Result, actual: *const Result) !void {
    switch (expected.*) {
        .atom => |expected_contacts| {
            const actual_contacts = actual.atom;
            try std.testing.expectEqual(expected_contacts.items.len, actual_contacts.items.len);
            for (expected_contacts.items, actual_contacts.items) |e, a| {
                try std.testing.expectEqual(e.a, a.a);
                try std.testing.expectEqual(e.b, a.b);
                try std.testing.expectEqual(e.distance, a.distance);
            }
        },
        .residue => |expected_contacts| {
            const actual_contacts = actual.residue;
            try std.testing.expectEqual(expected_contacts.items.len, actual_contacts.items.len);
            for (expected_contacts.items, actual_contacts.items) |e, a| {
                try std.testing.expectEqual(e.atom_a, a.atom_a);
                try std.testing.expectEqual(e.atom_b, a.atom_b);
                try std.testing.expectEqual(e.residue_a, a.residue_a);
                try std.testing.expectEqual(e.residue_b, a.residue_b);
                try std.testing.expectEqual(e.distance, a.distance);
            }
        },
    }
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

test "randomized cell list matches scalar reference across modes and selections" {
    const Case = struct { mode: Mode, cutoff: f64, select1: []const u8, select2: []const u8, scope: Scope };
    const cases = [_]Case{
        .{ .mode = .atom, .cutoff = 2.5, .select1 = "all", .select2 = "all", .scope = .all },
        .{ .mode = .atom, .cutoff = 4.0, .select1 = "heavy", .select2 = "heavy", .scope = .inter_residue },
        .{ .mode = .atom, .cutoff = 6.5, .select1 = "chain:A", .select2 = "chain:B", .scope = .inter_residue },
        .{ .mode = .atom, .cutoff = 8.0, .select1 = "protein,backbone", .select2 = "protein,sidechain", .scope = .all },
        .{ .mode = .residue, .cutoff = 2.5, .select1 = "all", .select2 = "all", .scope = .all },
        .{ .mode = .residue, .cutoff = 4.0, .select1 = "heavy", .select2 = "heavy", .scope = .inter_residue },
        .{ .mode = .residue, .cutoff = 6.5, .select1 = "chain:A", .select2 = "chain:B", .scope = .inter_residue },
        .{ .mode = .residue, .cutoff = 8.0, .select1 = "protein,backbone", .select2 = "protein,sidechain", .scope = .all },
    };
    const names = [_][]const u8{ "N", "CA", "CB", "O" };
    var prng = std.Random.DefaultPrng.init(0x5a434f4e54414354);
    const random = prng.random();
    for (0..24) |_| {
        var structure = model.Structure{};
        defer structure.deinit(std.testing.allocator);
        for (0..48) |i| {
            const residue_index: u32 = @intCast(i / 3);
            const atom_name = names[i % names.len];
            try structure.atoms.append(std.testing.allocator, .{
                .serial = @intCast(i + 1),
                .record = if (i % 11 == 0) .hetatm else .atom,
                .name = try model.Field.init(atom_name),
                .residue_name = try model.Field.init("ALA"),
                .chain = try model.Field.init(if (residue_index % 2 == 0) "A" else "B"),
                .residue_seq = @intCast(residue_index + 1),
                .element = try model.Field.init(if (i % 7 == 0) "H" else if (std.mem.eql(u8, atom_name, "N")) "N" else if (std.mem.eql(u8, atom_name, "O")) "O" else "C"),
                .x = random.float(f64) * 40.0 - 20.0,
                .y = random.float(f64) * 40.0 - 20.0,
                .z = random.float(f64) * 40.0 - 20.0,
                .residue_index = residue_index,
            });
        }
        for (cases) |case| {
            var expected = try calculateBruteForce(std.testing.allocator, &structure, case.mode, case.cutoff, case.select1, case.select2, case.scope);
            defer expected.deinit(std.testing.allocator);
            var actual = try calculate(std.testing.allocator, &structure, case.mode, case.cutoff, case.select1, case.select2, case.scope);
            defer actual.deinit(std.testing.allocator);
            try expectSameResult(&expected, &actual);
        }
    }
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
