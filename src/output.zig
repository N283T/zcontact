const std = @import("std");
const model = @import("model.zig");
const contact = @import("contact.zig");

pub fn writeResult(writer: *std.Io.Writer, structure: *const model.Structure, result: *const contact.Result) !void {
    switch (result.*) {
        .atom => |contacts| {
            try writer.writeAll("index1\tmodel1\tchain1\tresseq1\ticode1\tresname1\tatom1\taltloc1\tserial1\tindex2\tmodel2\tchain2\tresseq2\ticode2\tresname2\tatom2\taltloc2\tserial2\tdistance_A\n");
            for (contacts.items) |c| {
                const a = structure.atoms.items[c.a];
                const b = structure.atoms.items[c.b];
                try writer.print("{d}\t{d}\t{s}\t{d}\t{s}\t{s}\t{s}\t{s}\t{d}\t{d}\t{d}\t{s}\t{d}\t{s}\t{s}\t{s}\t{s}\t{d}\t{d:.3}\n", .{ c.a, a.model, a.chain.slice(), a.residue_seq, a.insertion.slice(), a.residue_name.slice(), a.name.slice(), a.altloc.slice(), a.serial, c.b, b.model, b.chain.slice(), b.residue_seq, b.insertion.slice(), b.residue_name.slice(), b.name.slice(), b.altloc.slice(), b.serial, c.distance });
            }
        },
        .residue => |contacts| {
            try writer.writeAll("residue_index1\tmodel1\tchain1\tresseq1\ticode1\tresname1\tresidue_index2\tmodel2\tchain2\tresseq2\ticode2\tresname2\tmin_distance_A\tatom1\taltloc1\tatom2\taltloc2\n");
            for (contacts.items) |c| {
                const a = structure.atoms.items[c.atom_a];
                const b = structure.atoms.items[c.atom_b];
                try writer.print("{d}\t{d}\t{s}\t{d}\t{s}\t{s}\t{d}\t{d}\t{s}\t{d}\t{s}\t{s}\t{d:.3}\t{s}\t{s}\t{s}\t{s}\n", .{ a.residue_index, a.model, a.chain.slice(), a.residue_seq, a.insertion.slice(), a.residue_name.slice(), b.residue_index, b.model, b.chain.slice(), b.residue_seq, b.insertion.slice(), b.residue_name.slice(), c.distance, a.name.slice(), a.altloc.slice(), b.name.slice(), b.altloc.slice() });
            }
        },
    }
}

pub fn writeSelectedAtoms(writer: *std.Io.Writer, structure: *const model.Structure, select1: []const u8, select2: []const u8) !void {
    const selection = @import("selection.zig");
    try writer.writeAll("index\tresidue_index\tmodel\tchain\tresseq\ticode\tresname\tatom\taltloc\telement\tx\ty\tz\n");
    for (structure.atoms.items, 0..) |atom, index| {
        if (!selection.matches(atom, select1) and !selection.matches(atom, select2)) continue;
        try writer.print("{d}\t{d}\t{d}\t{s}\t{d}\t{s}\t{s}\t{s}\t{s}\t{s}\t{d:.6}\t{d:.6}\t{d:.6}\n", .{ index, atom.residue_index, atom.model, atom.chain.slice(), atom.residue_seq, atom.insertion.slice(), atom.residue_name.slice(), atom.name.slice(), atom.altloc.slice(), atom.element.slice(), atom.x, atom.y, atom.z });
    }
}

pub fn count(result: *const contact.Result) usize {
    return switch (result.*) {
        .atom => |contacts| contacts.items.len,
        .residue => |contacts| contacts.items.len,
    };
}
