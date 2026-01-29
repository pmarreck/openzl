const std = @import("std");

const About = "Normally distributed random integer between bounds via Box-Muller sampling (deterministic seed optional)";

fn stdoutPrint(comptime fmt: []const u8, args: anytype) void {
    var buf: [1024]u8 = undefined;
    var w = std.fs.File.stdout().writer(&buf);
    _ = w.interface.print(fmt, args) catch {};
    _ = w.interface.flush() catch {};
}

fn stderrPrint(comptime fmt: []const u8, args: anytype) void {
    var buf: [1024]u8 = undefined;
    var w = std.fs.File.stderr().writer(&buf);
    _ = w.interface.print(fmt, args) catch {};
    _ = w.interface.flush() catch {};
}

const Options = struct {
    binary_output: bool = false,
    count: ?usize = null,
    start: ?i64 = null,
    end_val: ?i64 = null,
    seed: ?u64 = null,
};

fn showHelp() void {
    stdoutPrint(
        \\Usage: nrandom [options] [start] [end]
        \\Outputs normally-distributed random numbers between <start> and <end>
        \\If <start> is not specified, it defaults to 0
        \\If <end> is not specified, it defaults to 99
        \\
        \\Options:
        \\  -a, --about         Show a short description
        \\  -b, --binaryoutput  Output binary bytes (default 0-255, custom range allowed)
        \\  -c, --count N       Output N numbers (default: 1, or 1024 with -b)
        \\  -s, --seed N        Deterministic seed (decimal or hex, 0x prefix optional)
        \\  -h, --help          Show this help message
        \\
    , .{});
}

fn showAbout() void {
    stdoutPrint("{s}\n", .{About});
}

fn isHexDigit(c: u8) bool {
    return (c >= '0' and c <= '9') or (c >= 'a' and c <= 'f') or (c >= 'A' and c <= 'F');
}

fn parseSeed(s: []const u8) !u64 {
    if (std.mem.startsWith(u8, s, "0x") or std.mem.startsWith(u8, s, "0X")) {
        return std.fmt.parseUnsigned(u64, s[2..], 16);
    }
    var has_hex_letter = false;
    for (s) |ch| {
        if (!isHexDigit(ch)) {
            return std.fmt.parseUnsigned(u64, s, 10);
        }
        if (ch >= 'a' and ch <= 'f' or ch >= 'A' and ch <= 'F') {
            has_hex_letter = true;
        }
    }
    if (has_hex_letter) {
        stderrPrint(
            "Warning: seed contains hex letters but no 0x prefix; interpreting as hex\n",
            .{},
        );
        return std.fmt.parseUnsigned(u64, s, 16);
    }
    return std.fmt.parseUnsigned(u64, s, 10);
}

fn parseArgs() !Options {
    var opts = Options{};
    var pos1: ?[]const u8 = null;
    var pos2: ?[]const u8 = null;

    var it = try std.process.argsWithAllocator(std.heap.page_allocator);
    defer it.deinit();
    _ = it.next(); // skip argv[0]

    while (it.next()) |arg| {
        if (std.mem.eql(u8, arg, "--about") or std.mem.eql(u8, arg, "-a")) {
            showAbout();
            std.process.exit(0);
        } else if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            showHelp();
            std.process.exit(0);
        } else if (std.mem.eql(u8, arg, "--binaryoutput") or std.mem.eql(u8, arg, "-b")) {
            opts.binary_output = true;
        } else if (std.mem.eql(u8, arg, "--count") or std.mem.eql(u8, arg, "-c")) {
            const val = it.next() orelse {
                stderrPrint("Error: --count requires a number argument\n", .{});
                std.process.exit(1);
            };
            opts.count = std.fmt.parseUnsigned(usize, val, 10) catch {
                stderrPrint("Error: --count value must be a number\n", .{});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--seed") or std.mem.eql(u8, arg, "-s")) {
            const val = it.next() orelse {
                stderrPrint("Error: --seed requires a number argument\n", .{});
                std.process.exit(1);
            };
            opts.seed = parseSeed(val) catch {
                stderrPrint("Error: --seed value must be a number\n", .{});
                std.process.exit(1);
            };
        } else if (std.mem.startsWith(u8, arg, "-")) {
            stderrPrint("Error: unknown option {s}\n", .{arg});
            std.process.exit(1);
        } else {
            if (pos1 == null) {
                pos1 = arg;
            } else if (pos2 == null) {
                pos2 = arg;
            } else {
                stderrPrint("Error: too many positional arguments\n", .{});
                std.process.exit(1);
            }
        }
    }

    if (pos1) |p1| {
        opts.start = std.fmt.parseInt(i64, p1, 10) catch {
            stderrPrint("Error: start value must be a number\n", .{});
            std.process.exit(1);
        };
    }
    if (pos2) |p2| {
        opts.end_val = std.fmt.parseInt(i64, p2, 10) catch {
            stderrPrint("Error: end value must be a number\n", .{});
            std.process.exit(1);
        };
    }
    return opts;
}

const Pcg32 = struct {
    state: u64,
    inc: u64,

    pub fn init(seed: u64) Pcg32 {
        var rng = Pcg32{
            .state = 0,
            .inc = (seed << 1) | 1,
        };
        _ = rng.nextU32();
        rng.state += seed;
        _ = rng.nextU32();
        return rng;
    }

    pub fn nextU32(self: *Pcg32) u32 {
        const oldstate = self.state;
        self.state = oldstate * 6364136223846793005 + self.inc;
        const xorshifted = @as(u32, @truncate(((oldstate >> 18) ^ oldstate) >> 27));
        const rot: u5 = @truncate(oldstate >> 59);
        return (xorshifted >> rot) | (xorshifted << (@as(u5, 0) -% rot));
    }

    pub fn nextF64(self: *Pcg32) f64 {
        const u = self.nextU32();
        return (@as(f64, @floatFromInt(u)) + 1.0) / 4294967297.0;
    }
};

fn normalInt(rng: *Pcg32, start: i64, end_val: i64) i64 {
    if (start == end_val) return start;
    const range = end_val - start;
    while (true) {
        var u1v = rng.nextF64();
        if (u1v <= 0.0) u1v = rng.nextF64();
        const u2v = rng.nextF64();

        const z0 = std.math.sqrt(-2.0 * std.math.log(f64, std.math.e, u1v)) * std.math.cos(2.0 * std.math.pi * u2v);
        const rand_val = @as(f64, @floatFromInt(start))
            + z0 * (@as(f64, @floatFromInt(range)) / 6.0)
            + (@as(f64, @floatFromInt(range)) / 2.0);
        const rounded = std.math.floor(rand_val + 0.5);
        const result = @as(i64, @intFromFloat(rounded));
        if (result >= start and result <= end_val) return result;
    }
}

pub fn main() !void {
    var opts = try parseArgs();

    if (opts.binary_output) {
        opts.start = opts.start orelse 0;
        opts.end_val = opts.end_val orelse 255;
        opts.count = opts.count orelse 1024;
    } else {
        if (opts.start == null and opts.end_val == null) {
            stderrPrint("(with a start of 0 and an end of 99)\n", .{});
        }
        opts.start = opts.start orelse 0;
        opts.end_val = opts.end_val orelse 99;
        opts.count = opts.count orelse 1;
    }

    if (opts.start.? > opts.end_val.?) {
        stderrPrint("Error: start must be <= end\n", .{});
        std.process.exit(1);
    }

    if (opts.binary_output) {
        if (opts.start.? < 0 or opts.end_val.? > 255) {
            stderrPrint("Error: binary output range must be within 0..255\n", .{});
            std.process.exit(1);
        }
    }

    const now = std.time.nanoTimestamp();
    const seed = opts.seed orelse @as(u64, @truncate(@as(u128, @bitCast(now))));
    var rng = Pcg32.init(seed);

    var out_buf: [16384]u8 = undefined;
    var stdout = std.fs.File.stdout().writer(&out_buf);

    if (opts.binary_output) {
        var i: usize = 0;
        while (i < opts.count.?) : (i += 1) {
            const v = normalInt(&rng, opts.start.?, opts.end_val.?);
            try stdout.interface.writeByte(@as(u8, @intCast(v)));
        }
    } else {
        var i: usize = 0;
        while (i < opts.count.?) : (i += 1) {
            const v = normalInt(&rng, opts.start.?, opts.end_val.?);
            try stdout.interface.print("{d}\n", .{v});
        }
    }
    try stdout.interface.flush();
}
