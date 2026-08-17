const std = @import("std");
const zcontact = @import("zcontact");

test "PDB and mmCIF fixtures have equivalent selected contacts" {
    const allocator = std.testing.allocator;
    var pdb = try zcontact.parser.parse(allocator, @embedFile("mini.pdb"), .pdb, 1);
    defer pdb.deinit(allocator);
    var cif = try zcontact.parser.parse(allocator, @embedFile("mini.cif"), .mmcif, 1);
    defer cif.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 4), pdb.atoms.items.len);
    try std.testing.expectEqual(@as(usize, 3), cif.atoms.items.len);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), pdb.atoms.items[0].x, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), cif.atoms.items[0].x, 1e-12);

    var pdb_contacts = try zcontact.contact.calculate(allocator, &pdb, .atom, 4.0, "polymer,heavy", "polymer,heavy", .inter_residue);
    defer pdb_contacts.deinit(allocator);
    var cif_contacts = try zcontact.contact.calculate(allocator, &cif, .atom, 4.0, "polymer,heavy", "polymer,heavy", .inter_residue);
    defer cif_contacts.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 2), pdb_contacts.atom.items.len);
    try std.testing.expectEqual(pdb_contacts.atom.items.len, cif_contacts.atom.items.len);
    for (pdb_contacts.atom.items, cif_contacts.atom.items) |a, b|
        try std.testing.expectApproxEqAbs(a.distance, b.distance, 1e-12);
}
