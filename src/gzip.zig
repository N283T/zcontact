//! Bounded gzip input using Zig's standard-library DEFLATE implementation.

const std = @import("std");

pub const default_max_size: usize = 1024 * 1024 * 1024;

/// Reads a plain or gzip-compressed file. Gzip is detected by magic bytes, not
/// only by its extension. The caller owns the returned slice.
pub fn readFile(allocator: std.mem.Allocator, path: []const u8, max_size: usize) ![]u8 {
    const io = std.Io.Threaded.global_single_threaded.io();
    const data = try std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(max_size));
    errdefer allocator.free(data);
    if (data.len < 2 or data[0] != 0x1f or data[1] != 0x8b) return data;

    var source: std.Io.Reader = .fixed(data);
    var window: [std.compress.flate.max_window_len]u8 = undefined;
    var decompressor: std.compress.flate.Decompress = .init(&source, .gzip, &window);
    var output: std.Io.Writer.Allocating = .init(allocator);
    errdefer output.deinit();
    var limited = decompressor.reader.limited(.limited(max_size + 1), &.{});
    const count = try limited.interface.streamRemaining(&output.writer);
    if (count > max_size) return error.StreamTooLong;
    const result = try output.toOwnedSlice();
    allocator.free(data);
    return result;
}

test "gzip input is decompressed and bounded" {
    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 'H',  'e',  'l',  'l',  'o',
        ' ',  'w',  'o',  'r',  'l',  'd',  '\n', 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });
    const path = try tmp.dir.realPathFileAlloc(std.testing.io, "test.gz", std.testing.allocator);
    defer std.testing.allocator.free(path);
    const value = try readFile(std.testing.allocator, path, 64);
    defer std.testing.allocator.free(value);
    try std.testing.expectEqualStrings("Hello world\n", value);
    try std.testing.expectError(error.StreamTooLong, readFile(std.testing.allocator, path, 5));
}
