const std = @import("std");
const model = @import("model.zig");

pub const InputFormat = enum { pdb, mmcif };

pub fn detectFormat(path: []const u8, source: []const u8) !InputFormat {
    if (std.mem.endsWith(u8, path, ".pdb") or std.mem.endsWith(u8, path, ".ent") or std.mem.endsWith(u8, path, ".pdb.gz") or std.mem.endsWith(u8, path, ".ent.gz")) return .pdb;
    if (std.mem.endsWith(u8, path, ".cif") or std.mem.endsWith(u8, path, ".mmcif") or std.mem.endsWith(u8, path, ".cif.gz") or std.mem.endsWith(u8, path, ".mmcif.gz")) return .mmcif;
    const trimmed = std.mem.trimStart(u8, source, " \t\r\n");
    if (std.mem.startsWith(u8, trimmed, "data_") or std.mem.indexOf(u8, trimmed, "_atom_site.") != null) return .mmcif;
    if (std.mem.startsWith(u8, trimmed, "ATOM  ") or std.mem.startsWith(u8, trimmed, "HETATM") or std.mem.indexOf(u8, trimmed, "\nATOM  ") != null or std.mem.indexOf(u8, trimmed, "\nHETATM") != null) return .pdb;
    return error.UnknownInputFormat;
}

pub fn parse(allocator: std.mem.Allocator, source: []const u8, format: InputFormat, wanted_model: u32) !model.Structure {
    var structure = switch (format) {
        .pdb => try parsePdb(allocator, source, wanted_model),
        .mmcif => try parseMmcif(allocator, source, wanted_model),
    };
    errdefer structure.deinit(allocator);
    if (structure.atoms.items.len == 0) return error.NoAtomsForModel;
    try resolveAltlocs(allocator, &structure);
    try assignResidues(allocator, &structure);
    return structure;
}

fn field(line: []const u8, start: usize, end: usize) []const u8 {
    if (line.len <= start) return "";
    return std.mem.trim(u8, line[start..@min(end, line.len)], " ");
}

fn parseHybrid36(text: []const u8) !i64 {
    const trimmed = std.mem.trim(u8, text, " ");
    if (trimmed.len == 0) return error.InvalidHybrid36;
    if (trimmed[0] == '-' or std.ascii.isDigit(trimmed[0])) return std.fmt.parseInt(i64, trimmed, 10);
    const width: u32 = @intCast(text.len);
    var base_value: i64 = 0;
    const lowercase = std.ascii.isLower(trimmed[0]);
    for (trimmed) |character| {
        const digit: i64 = if (std.ascii.isDigit(character))
            character - '0'
        else if (std.ascii.isUpper(character))
            character - 'A' + 10
        else if (std.ascii.isLower(character))
            character - 'a' + 10
        else
            return error.InvalidHybrid36;
        base_value = base_value * 36 + digit;
    }
    const block: i64 = std.math.pow(i64, 36, width - 1);
    const decimal_limit: i64 = std.math.pow(i64, 10, width);
    return if (lowercase)
        base_value + 16 * block + decimal_limit
    else
        base_value - 10 * block + decimal_limit;
}

fn parsePdb(allocator: std.mem.Allocator, source: []const u8, wanted_model: u32) !model.Structure {
    var result = model.Structure{};
    errdefer result.deinit(allocator);
    var current_model: u32 = 1;
    var segment: u32 = 0;
    var has_models = false;
    var in_model = true;
    var lines = std.mem.splitScalar(u8, source, '\n');
    while (lines.next()) |raw| {
        const line = std.mem.trimEnd(u8, raw, "\r");
        if (std.mem.startsWith(u8, line, "MODEL ")) {
            has_models = true;
            in_model = true;
            current_model = try std.fmt.parseInt(u32, std.mem.trim(u8, field(line, 10, 14), " "), 10);
            segment = 0;
            continue;
        }
        if (std.mem.startsWith(u8, line, "ENDMDL")) {
            in_model = false;
            continue;
        }
        if (std.mem.startsWith(u8, line, "TER   ")) {
            segment += 1;
            continue;
        }
        const is_atom = std.mem.startsWith(u8, line, "ATOM  ");
        const is_het = std.mem.startsWith(u8, line, "HETATM");
        if ((!is_atom and !is_het) or (has_models and !in_model) or current_model != wanted_model or line.len < 54) continue;
        const raw_name = line[12..16];
        const name = try model.Field.init(field(line, 12, 16));
        const element_text = field(line, 76, 78);
        const chain_text = field(line, 21, 22);
        var internal_buffer: [32]u8 = undefined;
        const internal_text = try std.fmt.bufPrint(&internal_buffer, "{s}#{d}", .{ chain_text, segment });
        const x = try std.fmt.parseFloat(f64, field(line, 30, 38));
        const y = try std.fmt.parseFloat(f64, field(line, 38, 46));
        const z = try std.fmt.parseFloat(f64, field(line, 46, 54));
        if (!std.math.isFinite(x) or !std.math.isFinite(y) or !std.math.isFinite(z)) return error.NonFiniteCoordinate;
        const occupancy = if (field(line, 54, 60).len == 0) 1.0 else try std.fmt.parseFloat(f64, field(line, 54, 60));
        if (!std.math.isFinite(occupancy) or occupancy < 0) return error.InvalidOccupancy;
        const parsed_serial = try parseHybrid36(line[6..11]);
        if (parsed_serial < 0 or parsed_serial > std.math.maxInt(u32)) return error.InvalidAtomSerial;
        const parsed_residue = try parseHybrid36(line[22..26]);
        if (parsed_residue < std.math.minInt(i32) or parsed_residue > std.math.maxInt(i32)) return error.InvalidResidueSequence;
        try result.atoms.append(allocator, .{
            .serial = @intCast(parsed_serial),
            .model = current_model,
            .record = if (is_atom) .atom else .hetatm,
            .name = name,
            .altloc = try model.Field.init(field(line, 16, 17)),
            .residue_name = try model.Field.init(field(line, 17, 20)),
            .chain = try model.Field.init(chain_text),
            .internal_chain = try model.Field.init(internal_text),
            .residue_seq = @intCast(parsed_residue),
            .insertion = try model.Field.init(field(line, 26, 27)),
            .element = if (element_text.len > 0) try model.Field.init(element_text) else try model.inferPdbElement(raw_name, is_het),
            .x = x,
            .y = y,
            .z = z,
            .occupancy = occupancy,
        });
    }
    return result;
}

const Tokenizer = struct {
    source: []const u8,
    pos: usize = 0,

    fn next(self: *Tokenizer) ?[]const u8 {
        while (self.pos < self.source.len) {
            while (self.pos < self.source.len and std.ascii.isWhitespace(self.source[self.pos])) self.pos += 1;
            if (self.pos >= self.source.len) return null;
            if (self.source[self.pos] == '#') {
                while (self.pos < self.source.len and self.source[self.pos] != '\n') self.pos += 1;
                continue;
            }
            const start = self.pos;
            if (self.source[self.pos] == ';' and (self.pos == 0 or self.source[self.pos - 1] == '\n')) {
                self.pos += 1;
                if (self.pos < self.source.len and self.source[self.pos] == '\r') self.pos += 1;
                if (self.pos < self.source.len and self.source[self.pos] == '\n') self.pos += 1;
                const content = self.pos;
                while (self.pos < self.source.len) : (self.pos += 1) {
                    if (self.source[self.pos] == ';' and (self.pos == 0 or self.source[self.pos - 1] == '\n')) {
                        var end = self.pos;
                        if (end > content and self.source[end - 1] == '\n') end -= 1;
                        if (end > content and self.source[end - 1] == '\r') end -= 1;
                        self.pos += 1;
                        return self.source[content..end];
                    }
                }
                return self.source[content..self.pos];
            }
            if (self.source[self.pos] == '\'' or self.source[self.pos] == '"') {
                const quote = self.source[self.pos];
                self.pos += 1;
                const content = self.pos;
                while (self.pos < self.source.len) : (self.pos += 1) {
                    if (self.source[self.pos] == quote and (self.pos + 1 == self.source.len or std.ascii.isWhitespace(self.source[self.pos + 1]))) break;
                }
                const quoted_value = self.source[content..self.pos];
                if (self.pos < self.source.len) self.pos += 1;
                return quoted_value;
            }
            while (self.pos < self.source.len and !std.ascii.isWhitespace(self.source[self.pos])) self.pos += 1;
            return self.source[start..self.pos];
        }
        return null;
    }
};

fn findColumn(tags: []const []const u8, names: []const []const u8) ?usize {
    for (names) |name| for (tags, 0..) |tag, i| if (std.ascii.eqlIgnoreCase(tag, name)) return i;
    return null;
}

fn value(row: []const []const u8, col: ?usize, fallback: []const u8) []const u8 {
    return if (col) |i| if (i < row.len) row[i] else fallback else fallback;
}

fn preferredValue(row: []const []const u8, preferred: ?usize, fallback_col: ?usize, default: []const u8) []const u8 {
    const preferred_value = value(row, preferred, "");
    if (preferred_value.len != 0 and !std.mem.eql(u8, preferred_value, ".") and !std.mem.eql(u8, preferred_value, "?")) return preferred_value;
    return value(row, fallback_col, default);
}

fn parseMmcif(allocator: std.mem.Allocator, source: []const u8, wanted_model: u32) !model.Structure {
    var result = model.Structure{};
    errdefer result.deinit(allocator);
    var tok = Tokenizer{ .source = source };
    while (tok.next()) |token| {
        if (!std.ascii.eqlIgnoreCase(token, "loop_")) continue;
        var tags = std.ArrayListUnmanaged([]const u8).empty;
        defer tags.deinit(allocator);
        var first_value: ?[]const u8 = null;
        while (tok.next()) |item| {
            if (std.mem.startsWith(u8, item, "_")) try tags.append(allocator, item) else {
                first_value = item;
                break;
            }
        }
        if (tags.items.len == 0 or !std.mem.startsWith(u8, tags.items[0], "_atom_site.")) continue;
        const group_c = findColumn(tags.items, &.{"_atom_site.group_PDB"});
        const id_c = findColumn(tags.items, &.{"_atom_site.id"});
        const elem_c = findColumn(tags.items, &.{"_atom_site.type_symbol"});
        const auth_name_c = findColumn(tags.items, &.{"_atom_site.auth_atom_id"});
        const label_name_c = findColumn(tags.items, &.{"_atom_site.label_atom_id"});
        const alt_c = findColumn(tags.items, &.{"_atom_site.label_alt_id"});
        const auth_resn_c = findColumn(tags.items, &.{"_atom_site.auth_comp_id"});
        const label_resn_c = findColumn(tags.items, &.{"_atom_site.label_comp_id"});
        const auth_chain_c = findColumn(tags.items, &.{"_atom_site.auth_asym_id"});
        const label_chain_c = findColumn(tags.items, &.{"_atom_site.label_asym_id"});
        const auth_seq_c = findColumn(tags.items, &.{"_atom_site.auth_seq_id"});
        const label_seq_c = findColumn(tags.items, &.{"_atom_site.label_seq_id"});
        const ins_c = findColumn(tags.items, &.{"_atom_site.pdbx_PDB_ins_code"});
        const x_c = findColumn(tags.items, &.{"_atom_site.Cartn_x"});
        const y_c = findColumn(tags.items, &.{"_atom_site.Cartn_y"});
        const z_c = findColumn(tags.items, &.{"_atom_site.Cartn_z"});
        const occ_c = findColumn(tags.items, &.{"_atom_site.occupancy"});
        const model_c = findColumn(tags.items, &.{"_atom_site.pdbx_PDB_model_num"});
        if ((auth_name_c == null and label_name_c == null) or (auth_resn_c == null and label_resn_c == null) or
            (auth_chain_c == null and label_chain_c == null) or (auth_seq_c == null and label_seq_c == null) or
            x_c == null or y_c == null or z_c == null) return error.MissingAtomSiteColumn;

        var row = std.ArrayListUnmanaged([]const u8).empty;
        defer row.deinit(allocator);
        var item = first_value;
        while (item) |v| {
            if (std.mem.startsWith(u8, v, "_") or std.ascii.eqlIgnoreCase(v, "loop_") or std.ascii.eqlIgnoreCase(v, "stop_") or std.mem.startsWith(u8, v, "data_")) break;
            try row.append(allocator, v);
            if (row.items.len == tags.items.len) {
                const row_model = try std.fmt.parseInt(u32, value(row.items, model_c, "1"), 10);
                if (row_model == wanted_model) {
                    const atom_name = try model.Field.init(preferredValue(row.items, auth_name_c, label_name_c, ""));
                    const group = value(row.items, group_c, "ATOM");
                    const is_atom_group = std.ascii.eqlIgnoreCase(group, "ATOM");
                    if (!is_atom_group and !std.ascii.eqlIgnoreCase(group, "HETATM")) return error.InvalidAtomSiteGroup;
                    const label_chain = preferredValue(row.items, label_chain_c, auth_chain_c, "");
                    const x = try std.fmt.parseFloat(f64, value(row.items, x_c, ""));
                    const y = try std.fmt.parseFloat(f64, value(row.items, y_c, ""));
                    const z = try std.fmt.parseFloat(f64, value(row.items, z_c, ""));
                    if (!std.math.isFinite(x) or !std.math.isFinite(y) or !std.math.isFinite(z)) return error.NonFiniteCoordinate;
                    const occupancy_text = value(row.items, occ_c, "1");
                    const occupancy = if (std.mem.eql(u8, occupancy_text, ".") or std.mem.eql(u8, occupancy_text, "?")) 1.0 else try std.fmt.parseFloat(f64, occupancy_text);
                    if (!std.math.isFinite(occupancy) or occupancy < 0) return error.InvalidOccupancy;
                    try result.atoms.append(allocator, .{
                        .serial = std.fmt.parseInt(u32, value(row.items, id_c, "0"), 10) catch @intCast(result.atoms.items.len + 1),
                        .model = row_model,
                        .record = if (is_atom_group) .atom else .hetatm,
                        .name = atom_name,
                        .altloc = try model.Field.init(value(row.items, alt_c, "")),
                        .residue_name = try model.Field.init(preferredValue(row.items, auth_resn_c, label_resn_c, "")),
                        .chain = try model.Field.init(preferredValue(row.items, auth_chain_c, label_chain_c, "")),
                        .internal_chain = try model.Field.init(label_chain),
                        .residue_seq = try std.fmt.parseInt(i32, preferredValue(row.items, auth_seq_c, label_seq_c, ""), 10),
                        .insertion = try model.Field.init(value(row.items, ins_c, "")),
                        .element = if (elem_c != null) try model.Field.init(value(row.items, elem_c, "")) else try model.inferElement(atom_name.slice()),
                        .x = x,
                        .y = y,
                        .z = z,
                        .occupancy = occupancy,
                    });
                }
                row.clearRetainingCapacity();
            }
            item = tok.next();
        }
        if (row.items.len != 0) return error.IncompleteAtomSiteRow;
        return result;
    }
    return error.NoAtomSiteLoop;
}

fn altBetter(candidate: model.Atom, incumbent: model.Atom) bool {
    const ca = candidate.altloc.slice();
    const ia = incumbent.altloc.slice();
    if (ca.len == 0 and ia.len != 0) return true;
    if (ca.len != 0 and ia.len == 0) return false;
    if (candidate.occupancy != incumbent.occupancy) return candidate.occupancy > incumbent.occupancy;
    if (ca.len == 0) return false;
    if (std.mem.eql(u8, ca, "A") or std.mem.eql(u8, ia, "A")) return std.mem.eql(u8, ca, "A");
    return std.mem.lessThan(u8, ca, ia);
}

const AtomSiteKey = struct {
    model_num: u32,
    name: model.Field,
    internal_chain: model.Field,
    seq: i32,
    insertion: model.Field,
};

const ResidueKey = struct {
    model_num: u32,
    internal_chain: model.Field,
    seq: i32,
    insertion: model.Field,
};

const ConformerKey = struct { residue: ResidueKey, label: model.Field };
const ConformerChoice = struct { label: model.Field, score: f64 };

fn residueKey(atom: model.Atom) ResidueKey {
    return .{ .model_num = atom.model, .internal_chain = model.Atom.internalId(atom), .seq = atom.residue_seq, .insertion = atom.insertion };
}

fn labelPreferred(candidate: model.Field, incumbent: model.Field) bool {
    const c = candidate.slice();
    const i = incumbent.slice();
    if (std.mem.eql(u8, c, "A") or std.mem.eql(u8, i, "A")) return std.mem.eql(u8, c, "A");
    return std.mem.lessThan(u8, c, i);
}

fn resolveAltlocs(allocator: std.mem.Allocator, structure: *model.Structure) !void {
    // Choose one coherent non-blank conformer label per residue position using
    // summed occupancy evidence. Blank atoms are shared by every conformer.
    var scores = std.AutoHashMapUnmanaged(ConformerKey, f64).empty;
    defer scores.deinit(allocator);
    for (structure.atoms.items) |atom| {
        if (atom.altloc.len == 0) continue;
        const entry = try scores.getOrPut(allocator, .{ .residue = residueKey(atom), .label = atom.altloc });
        if (!entry.found_existing) entry.value_ptr.* = 0;
        entry.value_ptr.* += atom.occupancy;
    }
    var choices = std.AutoHashMapUnmanaged(ResidueKey, ConformerChoice).empty;
    defer choices.deinit(allocator);
    var score_iterator = scores.iterator();
    while (score_iterator.next()) |entry| {
        const choice = try choices.getOrPut(allocator, entry.key_ptr.residue);
        if (!choice.found_existing or entry.value_ptr.* > choice.value_ptr.score or
            (entry.value_ptr.* == choice.value_ptr.score and labelPreferred(entry.key_ptr.label, choice.value_ptr.label)))
            choice.value_ptr.* = .{ .label = entry.key_ptr.label, .score = entry.value_ptr.* };
    }

    const OrderedAtom = struct { atom: model.Atom, order: usize };
    var site_order = std.AutoHashMapUnmanaged(AtomSiteKey, usize).empty;
    defer site_order.deinit(allocator);
    for (structure.atoms.items, 0..) |atom, order| {
        const key = AtomSiteKey{ .model_num = atom.model, .name = atom.name, .internal_chain = model.Atom.internalId(atom), .seq = atom.residue_seq, .insertion = atom.insertion };
        const entry = try site_order.getOrPut(allocator, key);
        if (!entry.found_existing) entry.value_ptr.* = order;
    }
    var ordered = std.ArrayListUnmanaged(OrderedAtom).empty;
    defer ordered.deinit(allocator);
    var selected_sites = std.AutoHashMapUnmanaged(AtomSiteKey, usize).empty;
    defer selected_sites.deinit(allocator);
    for (structure.atoms.items) |atom| {
        if (atom.altloc.len != 0) {
            const choice = choices.get(residueKey(atom)) orelse continue;
            if (!model.Field.eql(atom.altloc, choice.label)) continue;
        }
        const key = AtomSiteKey{
            .model_num = atom.model,
            .name = atom.name,
            .internal_chain = model.Atom.internalId(atom),
            .seq = atom.residue_seq,
            .insertion = atom.insertion,
        };
        if (selected_sites.get(key)) |idx| {
            if (altBetter(atom, ordered.items[idx].atom)) ordered.items[idx].atom = atom;
        } else {
            try selected_sites.put(allocator, key, ordered.items.len);
            try ordered.append(allocator, .{ .atom = atom, .order = site_order.get(key).? });
        }
    }
    std.mem.sort(OrderedAtom, ordered.items, {}, struct {
        fn lessThan(_: void, a: OrderedAtom, b: OrderedAtom) bool {
            return a.order < b.order;
        }
    }.lessThan);
    var chosen = std.ArrayListUnmanaged(model.Atom).empty;
    errdefer chosen.deinit(allocator);
    try chosen.ensureTotalCapacity(allocator, ordered.items.len);
    for (ordered.items) |entry| chosen.appendAssumeCapacity(entry.atom);
    structure.atoms.deinit(allocator);
    structure.atoms = chosen;
}

fn assignResidues(allocator: std.mem.Allocator, structure: *model.Structure) !void {
    var residues = std.AutoHashMapUnmanaged(ResidueKey, u32).empty;
    defer residues.deinit(allocator);
    for (structure.atoms.items) |*atom| {
        const key = residueKey(atom.*);
        const entry = try residues.getOrPut(allocator, key);
        if (!entry.found_existing) entry.value_ptr.* = @intCast(residues.count() - 1);
        atom.residue_index = entry.value_ptr.*;
    }
}

test "PDB model and occupancy altloc selection" {
    const source =
        "MODEL        1\n" ++
        "ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.40 10.00           C  \n" ++
        "ATOM      2  CA BALA A   1       1.000   0.000   0.000  0.60 10.00           C  \n" ++
        "ENDMDL\nMODEL        2\n" ++
        "ATOM      3  CA  ALA A   1       9.000   0.000   0.000  1.00 10.00           C  \nENDMDL\n";
    var s = try parse(std.testing.allocator, source, .pdb, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), s.atoms.items.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1), s.atoms.items[0].x, 1e-9);
    var second = try parse(std.testing.allocator, source, .pdb, 2);
    defer second.deinit(std.testing.allocator);
    try std.testing.expectApproxEqAbs(@as(f64, 9), second.atoms.items[0].x, 1e-9);
}

test "non-contiguous alternate records resolve to their first site position" {
    const source =
        "ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.40 10.00           C  \n" ++
        "ATOM      2  CA  GLY A   2       5.000   0.000   0.000  1.00 10.00           C  \n" ++
        "ATOM      3  CA BALA A   1       1.000   0.000   0.000  0.60 10.00           C  \n";
    var s = try parse(std.testing.allocator, source, .pdb, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), s.atoms.items.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1), s.atoms.items[0].x, 1e-9);
    try std.testing.expectEqualStrings("B", s.atoms.items[0].altloc.slice());
}

test "alternate selection is residue-coherent rather than per-atom chimera" {
    const source =
        "ATOM      1  N  AALA A   1       1.000   0.000   0.000  0.40 10.00           N  \n" ++
        "ATOM      2  N  BALA A   1       2.000   0.000   0.000  0.60 10.00           N  \n" ++
        "ATOM      3  CB AALA A   1       3.000   0.000   0.000  0.90 10.00           C  \n" ++
        "ATOM      4  CB BALA A   1       4.000   0.000   0.000  0.10 10.00           C  \n";
    var s = try parse(std.testing.allocator, source, .pdb, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), s.atoms.items.len);
    try std.testing.expectEqualStrings("A", s.atoms.items[0].altloc.slice());
    try std.testing.expectEqualStrings("A", s.atoms.items[1].altloc.slice());
    try std.testing.expectApproxEqAbs(@as(f64, 1), s.atoms.items[0].x, 1e-9);
}

test "duplicate blank atom sites keep highest occupancy" {
    const source =
        "ATOM      1  CA  ALA A   1       1.000   0.000   0.000  0.60 10.00           C  \n" ++
        "ATOM      2  CA  ALA A   1       0.000   0.000   0.000  0.40 10.00           C  \n";
    var s = try parse(std.testing.allocator, source, .pdb, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 1), s.atoms.items.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1), s.atoms.items[0].x, 1e-9);
}

test "TER separates repeated blank-chain residue identifiers" {
    const source =
        "ATOM      1  CA  ALA     1       0.000   0.000   0.000  1.00 10.00           C  \n" ++
        "TER       2      ALA     1                                                      \n" ++
        "ATOM      3  CA  GLY     1       3.000   0.000   0.000  1.00 10.00           C  \n";
    var s = try parse(std.testing.allocator, source, .pdb, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), s.atoms.items.len);
    try std.testing.expect(s.atoms.items[0].residue_index != s.atoms.items[1].residue_index);
}

test "mmCIF prefers author identifiers and keeps insertion codes distinct" {
    const source =
        "data_test\nloop_\n" ++
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" ++
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n" ++
        "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n" ++
        "_atom_site.auth_atom_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n" ++
        "ATOM 1 C CA ALA X 1 A 0 0 0 CA ALA A 10\n" ++
        "ATOM 2 C CA GLY X 2 B 2 0 0 CA GLY A 10\n#\n";
    var s = try parse(std.testing.allocator, source, .mmcif, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), s.atoms.items.len);
    try std.testing.expectEqualStrings("A", s.atoms.items[0].chain.slice());
    try std.testing.expectEqual(@as(i32, 10), s.atoms.items[0].residue_seq);
    try std.testing.expect(s.atoms.items[0].residue_index != s.atoms.items[1].residue_index);
}

test "mmCIF falls back from null author identifiers row by row" {
    const source =
        "data_test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n" ++
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n" ++
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n" ++
        "_atom_site.auth_atom_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n" ++
        "HETATM 1 O O HOH W 7 0 0 0 ? ? ? ?\n#\n";
    var s = try parse(std.testing.allocator, source, .mmcif, 1);
    defer s.deinit(std.testing.allocator);
    try std.testing.expectEqualStrings("W", s.atoms.items[0].chain.slice());
    try std.testing.expectEqual(@as(i32, 7), s.atoms.items[0].residue_seq);
    try std.testing.expectEqualStrings("HOH", s.atoms.items[0].residue_name.slice());
}

test "gzip-style suffix and HETATM content detection" {
    try std.testing.expectEqual(InputFormat.pdb, try detectFormat("ligand.pdb.gz", ""));
    try std.testing.expectEqual(InputFormat.mmcif, try detectFormat("structure.cif.gz", ""));
    try std.testing.expectEqual(InputFormat.pdb, try detectFormat("unknown", "HETATM    1"));
}

test "PDB hybrid-36 identifiers" {
    try std.testing.expectEqual(@as(i64, 9999), try parseHybrid36("9999"));
    try std.testing.expectEqual(@as(i64, 10000), try parseHybrid36("A000"));
    try std.testing.expectEqual(@as(i64, 1223055), try parseHybrid36("ZZZZ"));
    try std.testing.expectEqual(@as(i64, 1223056), try parseHybrid36("a000"));
}
