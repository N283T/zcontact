const std = @import("std");
const model = @import("model.zig");

pub fn matches(atom: model.Atom, expression: []const u8) bool {
    var clauses = std.mem.splitScalar(u8, expression, ',');
    while (clauses.next()) |raw| {
        const clause = std.mem.trim(u8, raw, " ");
        if (clause.len == 0 or std.ascii.eqlIgnoreCase(clause, "all")) continue;
        if (std.ascii.eqlIgnoreCase(clause, "polymer")) {
            if (atom.record != .atom) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "protein")) {
            if (atom.record != .atom or !isAminoAcid(atom.residue_name.slice())) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "hetero")) {
            if (atom.record != .hetatm) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "water")) {
            if (!isWater(atom.residue_name.slice())) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "ligand")) {
            if (atom.record != .hetatm or isWater(atom.residue_name.slice())) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "heavy")) {
            if (std.ascii.eqlIgnoreCase(atom.element.slice(), "H") or std.ascii.eqlIgnoreCase(atom.element.slice(), "D")) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "backbone")) {
            const n = atom.name.slice();
            if (!isBackbone(n)) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "sidechain")) {
            if (isBackbone(atom.name.slice())) return false;
        } else if (std.ascii.eqlIgnoreCase(clause, "ca")) {
            if (!std.mem.eql(u8, atom.name.slice(), "CA")) return false;
        } else if (splitPredicate(clause, "chain:")) |v| {
            const wanted = if (std.mem.eql(u8, v, "_")) "" else v;
            if (!std.mem.eql(u8, atom.chain.slice(), wanted)) return false;
        } else if (splitPredicate(clause, "name:")) |v| {
            if (!std.ascii.eqlIgnoreCase(atom.name.slice(), v)) return false;
        } else if (splitPredicate(clause, "element:")) |v| {
            if (!std.ascii.eqlIgnoreCase(atom.element.slice(), v)) return false;
        } else if (splitPredicate(clause, "resname:")) |v| {
            if (!std.ascii.eqlIgnoreCase(atom.residue_name.slice(), v)) return false;
        } else if (splitPredicate(clause, "resseq:")) |v| {
            const wanted = std.fmt.parseInt(i32, v, 10) catch return false;
            if (atom.residue_seq != wanted) return false;
        } else if (splitPredicate(clause, "icode:")) |v| {
            const wanted = if (std.mem.eql(u8, v, "_")) "" else v;
            if (!std.mem.eql(u8, atom.insertion.slice(), wanted)) return false;
        } else return false;
    }
    return true;
}

pub fn validate(expression: []const u8) !void {
    if (std.mem.trim(u8, expression, " ").len == 0 or std.mem.indexOfAny(u8, expression, "\t\r\n") != null) return error.InvalidSelection;
    var clauses = std.mem.splitScalar(u8, expression, ',');
    while (clauses.next()) |raw| {
        const c = std.mem.trim(u8, raw, " ");
        if (c.len == 0) return error.InvalidSelection;
        if (std.ascii.eqlIgnoreCase(c, "all") or std.ascii.eqlIgnoreCase(c, "polymer") or std.ascii.eqlIgnoreCase(c, "protein") or
            std.ascii.eqlIgnoreCase(c, "heavy") or std.ascii.eqlIgnoreCase(c, "backbone") or std.ascii.eqlIgnoreCase(c, "sidechain") or
            std.ascii.eqlIgnoreCase(c, "ca") or std.ascii.eqlIgnoreCase(c, "water") or std.ascii.eqlIgnoreCase(c, "hetero") or std.ascii.eqlIgnoreCase(c, "ligand")) continue;
        if (splitPredicate(c, "chain:") != null or splitPredicate(c, "name:") != null or
            splitPredicate(c, "element:") != null or splitPredicate(c, "resname:") != null or splitPredicate(c, "icode:") != null) continue;
        if (splitPredicate(c, "resseq:")) |v| {
            _ = std.fmt.parseInt(i32, v, 10) catch return error.InvalidSelection;
            continue;
        }
        return error.InvalidSelection;
    }
}

fn isAminoAcid(name: []const u8) bool {
    const names = [_][]const u8{
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "SEC", "PYL",
    };
    for (names) |candidate| if (std.ascii.eqlIgnoreCase(name, candidate)) return true;
    return false;
}

fn isWater(name: []const u8) bool {
    return std.ascii.eqlIgnoreCase(name, "HOH") or std.ascii.eqlIgnoreCase(name, "WAT") or
        std.ascii.eqlIgnoreCase(name, "H2O") or std.ascii.eqlIgnoreCase(name, "DOD");
}

fn isBackbone(name: []const u8) bool {
    return std.mem.eql(u8, name, "N") or std.mem.eql(u8, name, "CA") or std.mem.eql(u8, name, "C") or std.mem.eql(u8, name, "O");
}

fn splitPredicate(clause: []const u8, prefix: []const u8) ?[]const u8 {
    if (clause.len > prefix.len and std.ascii.startsWithIgnoreCase(clause, prefix)) return clause[prefix.len..];
    return null;
}

test "selection conjunction" {
    const atom = model.Atom{ .serial = 1, .name = try model.Field.init("CA"), .residue_name = try model.Field.init("ALA"), .chain = try model.Field.init("A"), .residue_seq = 1, .element = try model.Field.init("C"), .x = 0, .y = 0, .z = 0 };
    try std.testing.expect(matches(atom, "polymer,protein,heavy,backbone,chain:A"));
    try std.testing.expect(!matches(atom, "element:N"));
    var blank = atom;
    blank.chain = try model.Field.init("");
    try std.testing.expect(matches(blank, "chain:_"));
    try std.testing.expect(matches(atom, "resseq:1,icode:_"));
    var water = atom;
    water.record = .hetatm;
    water.residue_name = try model.Field.init("HOH");
    try std.testing.expect(matches(water, "water,hetero"));
    try std.testing.expect(!matches(water, "ligand"));
    water.residue_name = try model.Field.init("ATP");
    try std.testing.expect(matches(water, "ligand,heavy"));
}
