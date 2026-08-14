#!/usr/bin/env perl
#
# Generates the golden vectors that pin `plethora-compat`'s RANDLIB port to the
# Perl Math::Random that plethora itself calls.
#
# The chain reproduced here is exactly the one merge_pairs.pl walks:
#
#     md5_hex($line)  ->  random_set_seed_from_phrase($eed)  ->  random_normal(1, $mean, $sd)
#                          (= salfph = phrtsd + setall)          (= gennor = sd * snorm() + av)
#
# Run through the project-local oracle library so the system perl5 tree, which
# holds XS modules built against another perl, cannot interfere:
#
#     env -u PERL5LIB perl -I.oracle/perl5/lib/perl5 \
#         crates/plethora-compat/tests/oracle/gen_randlib_vectors.pl \
#         > crates/plethora-compat/tests/data/randlib_vectors.tsv
#
# %.17g is the shortest round-trip form for an IEEE double, so the Rust side can
# compare bit patterns rather than tolerances.

use strict;
use warnings;
use Math::Random qw(random_seed_from_phrase random_set_seed_from_phrase
                    random_uniform random_normal);
use Digest::MD5 qw(md5_hex);

# Phrases split three ways: what plethora really passes (32 hex characters from
# md5_hex), the edges of phrtsd's own character table, and the degenerate cases.
my @bedpe_lines = (
    "chr1\t4609301\t4609401\tchr1\t4609650\t4609750\tread1\t42\t+\t-",
    "chr1\t145289900\t145290000\tchr1\t145290400\t145290500\tsim_read_00017\t255\t-\t+",
    ".\t-1\t-1\tchr1\t4610000\t4610100\torphan\t0\t.\t+",
);

my @phrases = map { md5_hex($_) } @bedpe_lines;
push @phrases,
    q{},                                   # empty: phrtsd returns the defaults untouched
    q{ },                                  # a lone blank: lennob trims it to length 0
    q{a},                                  # one character
    q{ab},                                 # two: exercises the new phrtsd's dropped last char
    q{abc},
    q{trailing blanks   },                 # lennob strips these, interior blanks stay
    q{!@#$%^&*()_+[];:'"<>?,./},           # every punctuation glyph in the table
    q{ABCDEFGHIJKLMNOPQRSTUVWXYZ},
    q{0123456789},
    q{~|}. chr(1),                         # characters absent from the table
    "\x{c3}\x{a9}",                        # a UTF-8 sequence arriving as two bytes
    q{d41d8cd98f00b204e9800998ecf8427e},   # md5 of the empty string
    ('x' x 200);                           # longer than any md5

print "# phrase_hex\tseed1\tseed2\tuniform_x5\tnormal01_x5\tgennor_317_45\n";

for my $phrase (@phrases) {
    # Hex-encode the phrase so the TSV stays one record per line whatever the
    # phrase contains (tabs, quotes, high bytes).
    my $phrase_hex = unpack 'H*', $phrase;

    # phrtsd alone, no generator state touched.
    my @seed = random_seed_from_phrase($phrase);

    # ranf() straight through: random_uniform's default range is (0, 1) and
    # genunf(low, high) is low + (high - low) * ranf().
    random_set_seed_from_phrase($phrase);
    my @uniform = map { random_uniform() } 1 .. 5;

    # snorm() straight through: gennor(0, 1) is 1 * snorm() + 0.
    random_set_seed_from_phrase($phrase);
    my @normal = map { random_normal() } 1 .. 5;

    # The shape merge_pairs.pl uses: one draw, mean and sd from the sample.
    random_set_seed_from_phrase($phrase);
    my $gennor = random_normal(1, 317, 45);

    print join("\t",
        $phrase_hex,
        $seed[0],
        $seed[1],
        join(',', map { sprintf '%.17g', $_ } @uniform),
        join(',', map { sprintf '%.17g', $_ } @normal),
        sprintf('%.17g', $gennor),
    ), "\n";
}
