const std = @import("std");

pub const Field = struct {
    bytes: [32]u8 = [_]u8{0} ** 32,
    len: u8 = 0,

    pub fn init(value: []const u8) !Field {
        const v = if (std.mem.eql(u8, value, ".") or std.mem.eql(u8, value, "?")) "" else value;
        if (v.len > 32) return error.FieldTooLong;
        if (std.mem.indexOfAny(u8, v, "\t\r\n") != null) return error.InvalidFieldCharacter;
        var result = Field{};
        @memcpy(result.bytes[0..v.len], v);
        result.len = @intCast(v.len);
        return result;
    }

    pub fn slice(self: *const Field) []const u8 {
        return self.bytes[0..self.len];
    }

    pub fn eql(a: Field, b: Field) bool {
        return std.mem.eql(u8, a.slice(), b.slice());
    }
};

pub const Atom = struct {
    serial: u32,
    model: u32 = 1,
    record: enum { atom, hetatm } = .atom,
    name: Field,
    altloc: Field = .{},
    residue_name: Field,
    chain: Field,
    /// Parser-internal molecular/segment identity. PDB increments this at TER;
    /// mmCIF uses label_asym_id. Output and selections continue to use `chain`.
    internal_chain: Field = .{},
    residue_seq: i32,
    insertion: Field = .{},
    element: Field,
    x: f64,
    y: f64,
    z: f64,
    occupancy: f64 = 1.0,
    residue_index: u32 = 0,

    pub fn sameSite(a: Atom, b: Atom) bool {
        return a.model == b.model and a.residue_seq == b.residue_seq and
            Field.eql(internalId(a), internalId(b)) and Field.eql(a.insertion, b.insertion) and
            Field.eql(a.residue_name, b.residue_name) and Field.eql(a.name, b.name);
    }

    pub fn sameResidue(a: Atom, b: Atom) bool {
        return a.model == b.model and a.residue_seq == b.residue_seq and
            Field.eql(internalId(a), internalId(b)) and Field.eql(a.insertion, b.insertion);
    }

    pub fn internalId(atom: Atom) Field {
        return if (atom.internal_chain.len != 0) atom.internal_chain else atom.chain;
    }
};

pub const Structure = struct {
    atoms: std.ArrayListUnmanaged(Atom) = .empty,

    pub fn deinit(self: *Structure, allocator: std.mem.Allocator) void {
        self.atoms.deinit(allocator);
    }
};

pub fn inferElement(atom_name: []const u8) !Field {
    var name = std.mem.trim(u8, atom_name, " ");
    while (name.len > 0 and std.ascii.isDigit(name[0])) name = name[1..];
    if (name.len == 0) return Field.init("");
    var one: [1]u8 = .{std.ascii.toUpper(name[0])};
    return Field.init(&one);
}

pub fn inferPdbElement(raw_name: []const u8, is_hetatm: bool) !Field {
    if (raw_name.len < 2) return inferElement(raw_name);
    if (raw_name[0] == ' ' or std.ascii.isDigit(raw_name[0])) {
        var one: [1]u8 = .{std.ascii.toUpper(raw_name[1])};
        return Field.init(&one);
    }
    if (is_hetatm and std.ascii.isAlphabetic(raw_name[1])) {
        var two: [2]u8 = .{ std.ascii.toUpper(raw_name[0]), std.ascii.toUpper(raw_name[1]) };
        return Field.init(&two);
    }
    var one: [1]u8 = .{std.ascii.toUpper(raw_name[0])};
    return Field.init(&one);
}

test "PDB element inference respects atom-name alignment" {
    try std.testing.expectEqualStrings("C", (try inferPdbElement(" CA ", false)).slice());
    try std.testing.expectEqualStrings("FE", (try inferPdbElement("FE  ", true)).slice());
    try std.testing.expectEqualStrings("CL", (try inferPdbElement("CL  ", true)).slice());
    try std.testing.expectEqualStrings("H", (try inferPdbElement("1HG ", false)).slice());
}
