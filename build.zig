const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("zcontact", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });
    const options = b.addOptions();
    options.addOption([]const u8, "version", "0.1.0");
    const exe = b.addExecutable(.{
        .name = "zcontact",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zcontact", .module = mod },
                .{ .name = "build_options", .module = options.createModule() },
            },
        }),
    });
    b.installArtifact(exe);

    const run_step = b.step("run", "Run zcontact");
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| run_cmd.addArgs(args);
    run_step.dependOn(&run_cmd.step);

    const tests = b.addTest(.{ .root_module = mod });
    const integration_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/integration.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zcontact", .module = mod }},
        }),
    });
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&b.addRunArtifact(tests).step);
    test_step.dependOn(&b.addRunArtifact(integration_tests).step);
    const cli_tests = b.addSystemCommand(&.{"sh"});
    cli_tests.addFileArg(b.path("tests/cli-smoke.sh"));
    cli_tests.addArtifactArg(exe);
    cli_tests.addFileArg(b.path("tests/mini.pdb"));
    test_step.dependOn(&cli_tests.step);

    const bench = b.addExecutable(.{
        .name = "zcontact-bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/benchmark.zig"),
            .target = target,
            .optimize = .ReleaseFast,
            .imports = &.{.{ .name = "zcontact", .module = mod }},
        }),
    });
    const bench_cmd = b.addRunArtifact(bench);
    const bench_step = b.step("bench", "Run deterministic cell-list microbenchmarks");
    bench_step.dependOn(&bench_cmd.step);
}
