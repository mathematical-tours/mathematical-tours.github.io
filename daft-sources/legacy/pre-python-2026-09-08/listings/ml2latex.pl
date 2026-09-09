#!/usr/bin/env perl -w
#
# ml2latex - Matlab to LaTeX converter.
# Copyright (C) 1998  Peter J. Acklam
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# The GNU General Public License can be found on the World Wide Web at
# http://www.gnu.org/copyleft/gpl.html and it can be obtained by writing
# to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#
# Author:      Peter J. Acklam
# Time-stamp:  2001-12-15 23:53:08 +0100
# E-mail:      pjacklam@online.no
# URL:         http://home.online.no/~pjacklam
#
# Description
#
# This is a Matlab to LaTeX converter. It takes names of Matlab files
# from the command line or the standard input and generates the
# corresponding LaTeX file with the same name, but a .tex suffix. It can
# also convert Matlab code from the standard input and write the
# corresponding LaTeX code on the standard output.
#
# Installation
#
# Save this file to disk and make it executable. Edit the top line so it
# points to your local Perl interpreter. Use the -h or --help switch to
# display program usage.
#
# Things to do...perhaps...some day...
#
#   o Write a manual page.
#   o Make program perform filename globbing.
#   o Make program recurse subdirectories.
#   o Make program write index of Matlab subfunctions.
#   o Get more lines on each page.
#   o Allow different font sizes.
#   o Use delayed writing to avoid trailing empty lines in TeX file.
#   o Better handling of implicit comments (anything after ... ).

########################################################################
# Miscellaneous "high-level" stuff.
########################################################################

require 5.004;                  # Required version of Perl.

# Pragmatic modules.
use strict;                     # Restrict unsafe constructs.
#use diagnostics;                # Force verbose warning diagnostics.
use vars qw( $VERSION );        # Predeclare global variable names.

# Modules from the Standard Perl Library.
use File::Basename;     # Split a pathname into pieces.
use File::Path;         # Create or remove a series of directories.
#use File::Find;         # Traverse a file tree.
#use File::Copy;         # Copy files or filehandles.
use File::Spec;         # Portably perform operations on file names.
use Text::Tabs;         # Expand and unexpand tabs.
use Getopt::Long;       # Extended processing of command line options.
use FileHandle;         # Object methods for filehandles.
#use DirHandle;          # Object methods for directory handles.
#use Cwd;                # Get pathname of current working directory.

# Global variables.
$VERSION = "0.010";

# Lexical (private) variables.
my $program = fileparse $0;      # Name of this program.
my $is_ms   = $^O =~ /^(MS)?(Win32|DOS)$/i;

my $ml_normal_regex;    # Regular expression defined below.
my $ml_string_regex;    # Regular expression defined below.
my $ml_comment_regex;   # Regular expression defined below.
my $ml_keyword_regex;   # Regular expression defined below.

# Some default values.
my $begin_normal  = '\texttt{';
my $begin_string  = '\texttt{\textit{';
my $begin_comment = '\textit{';
my $begin_keyword = '\texttt{\textbf{';

########################################################################
# print_version
#
# Print program version and copyright information.
#
sub print_version () {
    print <<"EOF";

$program $VERSION

Copyright (c) 1998-1999 Peter J. Acklam. All rights reserved.
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License.

EOF
}

########################################################################
# print_usage
#
# Print program usage.
#
sub print_usage () {
    print <<"EOF";
Usage: $program [-fhlnopsvV] [-d DIR] [files...]
Matlab to LaTeX converter.

  -d, --outputdir DIR       output to given directory, not m-file dir
  -f, --fast                skip initialization file ~/.${program}rc
  -h, --help                print usage and exit
  -l, --linenum             include line numbers
  -n, --names               read names (not code) from standard input
  -o, --online              do the online help only
  -p, --pretty              include heading with date, filename and page
                            number
  -s, --standalone          generate stand-alone LaTeX document
  -v, --verbose             verbose mode
  -V, --version             print version information and exit
      --begin-normal STR    LaTeX code preceding "normal" Matlab code
                            (default is `$begin_normal')
      --begin-string STR    LaTeX code preceding Matlab strings, default
                            is `$begin_string'
      --begin-comment STR   LaTeX code preceding Matlab comments,
                            default is `$begin_comment'
      --begin-keyword STR   LaTeX code preceding Matlab keywords,
                            default is `$begin_keyword'

Report bugs to jacklam\@math.uio.no
EOF
}

########################################################################
# print_manual
#
# Print the manual page (POD) embedded in this file.
#
sub print_manual () {
    my $cmd = "perldoc $0";
    system( $cmd ) == 0 or die "$program: system ``$cmd'' failed\n";
}

########################################################################
# Define option variables.  They are private variables global to this
# file.  Configure option handling and get the user options.
#
#Getopt::Long::Configure( "no_ignore_case", "no_auto_abbrev", "bundling" );
Getopt::Long::config( "no_ignore_case", "no_auto_abbrev", "bundling" );
GetOptions( "d|outputdir=s"   => \my $opt_outputdir,
            "f|fast"          => \my $opt_fast,
            "l|linenum"       => \my $opt_linenum,
            "n|names"         => \my $opt_names,
            "o|online"        => \my $opt_online,
            "p|pretty"        => \my $opt_pretty,
            "s|standalone"    => \my $opt_standalone,
            "begin-normal=s"  => \my $opt_begin_normal,
            "begin-string=s"  => \my $opt_begin_string,
            "begin-comment=s" => \my $opt_begin_comment,
            "begin-keyword=s" => \my $opt_begin_keyword,

            "h|help"           => \my $opt_help,
            "m|manual"         => \my $opt_manual,
            "v|verbose"        => \my $opt_verbose,
            "V|version"        => \my $opt_version,
  )
  or die <<"EOF";
$program: usage error
Try `$program --help' for more information.
EOF

########################################################################
# Some options result in a quick exit.
#
print_usage,   exit if $opt_help;
print_manual,  exit if $opt_manual;
print_version, exit if $opt_version;

# Assign default values to options not given.
$begin_normal  = defined $opt_begin_normal  ? $opt_begin_normal  : '\texttt{';
$begin_string  = defined $opt_begin_string  ? $opt_begin_string  : '\texttt{\textit{';
$begin_comment = defined $opt_begin_comment ? $opt_begin_comment : '\textit{';
$begin_keyword = defined $opt_begin_keyword ? $opt_begin_keyword : '\texttt{\textbf{';

# Use left brace rather than empty string.
$begin_normal  = '{' if $begin_normal  eq '';
$begin_string  = '{' if $begin_string  eq '';
$begin_comment = '{' if $begin_comment eq '';
$begin_keyword = '{' if $begin_keyword eq '';

# The corresponding closing string contains only right braces. Make sure
# the number of right braces in the closing string matches the number of
# left braces in the opening string.
my $end_normal  = '}' x $begin_normal  =~ tr/{/{/;
my $end_string  = '}' x $begin_string  =~ tr/{/{/;
my $end_comment = '}' x $begin_comment =~ tr/{/{/;
my $end_keyword = '}' x $begin_keyword =~ tr/{/{/;

# The output directory.
my $outputdir = $opt_outputdir;
$opt_outputdir = defined $opt_outputdir ? 1 : 0;
mkpath $outputdir or die "$program: can't create $outputdir: $!\n"
  if $opt_outputdir && ! -d $outputdir;

########################################################################
# read_rcfile
#
# Read the initialization file and prepend data to argument list.
#
sub read_rcfile () {

    return if $opt_fast;

    # Get the home directory.
    my $home = $ENV{HOME} || $ENV{LOGDIR}
      || $^O !~ /^MS/ && (getpwuid($<))[7]
        || do {
               warn <<"EOF";
$program: home dir not found (HOME not set) so skipping ~/.${program}rc
EOF
               return;
           };

    # Find initialization file.
    my $rcfile = File::Spec->catfile( $home, ".${program}rc" );
    # Allow init file name to begin with "_" rather than "." when under
    # certain operating systems.
    if ( $is_ms & ! -f $rcfile ) {
        $rcfile = File::Spec->catfile( $home, "_${program}rc" )
    }
    return unless -f $rcfile;

    # Read initialization file.
    my @rcdata = ();
    my $fh = FileHandle->new($rcfile)
      or die "$program: can't open $rcfile: $!\n";
    while ( <$fh> ) {
        next if /^\s*[!#]/;     # Skip lines starting with ! or #.
        push @rcdata, split;
    }
    $fh->close;

    # Prepend the initialization file data to the argument list.
    unshift @ARGV, @rcdata;

}

########################################################################
# transl_tex_char
#
# Ttranslate TeX special characters to literal characters.
#
sub transl_tex_char (@) {

    # Get input.
    my @chunk = @_ ? @_ : $_;

    # Translate characters.
    foreach my $chunk ( @chunk ) {

        $chunk = expand $chunk;                 # Expand tabs to spaces.
        my @parts = split /(\s+)/, $chunk;
        foreach my $part ( @parts ) {

            #
            # Fix spaces. Neither \space nor \obeyspaces seems to work
            # with T1 encoding, but \phantom seems to always work.
            #
            if ( $part =~ /^\s+$/ ) {
                $part = "\\phantom{" . ( "x" x length $part ) . "}";

            #
            # Split on non-word characters including underscore. The ten
            # TeX special characters #$%&\^_{}~ require special
            # treatment. Also, put braces around all other non-word
            # characters to avoid ligatures.
            #
            } else {
                my @sub_parts = split /([\W_])/, $part;
                foreach ( @sub_parts ) {
                    next if $_ eq "";
                    next if s/ ([\$#%&{}_]) /\\$1\{}/x;
                    next if s/ \\ /\\textbackslash{}/x;
                    next if s/ \^ /\\textasciicircum{}/x;
                    next if s/  ~ /\\textasciitilde{}/x;
                    s/(\W)/{$1}/;
                }
                $part = join "", @sub_parts;

            }   # of if

        }   # of foreach
        $chunk = join "", @parts;

        #
        # Convert \foo{}\bar to \foo\bar etc. This is really not
        # necessary but it reduces the size of the output file. The {}
        # after any \^ and \~ must not be removed, however.
        #
        $chunk =~ s/ ( [^\\^~] ) {} \\ /$1\\/gx;

    }

    # Return output.
    return          @chunk if         wantarray;    # List context.
    return join "", @chunk if defined wantarray;    # Scalar context.
    $_ =   join "", @chunk;                         # Void context.

}

########################################################################
# write_preamble
#
# Write the LaTeX document preamble.
#
sub write_preamble (*$) {

    my ( $handle, $filename ) = @_;

    if ( $opt_standalone ) {
        print $handle <<"EOF";
\\documentclass[a4paper,10pt,twoside]{article}
\\usepackage[latin1]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{times}
EOF
        if ( $opt_pretty ) {
            print $handle <<"EOF";
\\usepackage{fancyhdr}
EOF
        }
        print $handle <<"EOF";

\\setlength{\\oddsidemargin}{-1in}
\\addtolength{\\oddsidemargin}{2cm}
\\setlength{\\evensidemargin}{\\oddsidemargin}
\\setlength{\\textwidth}{\\paperwidth}
\\addtolength{\\textwidth}{-2in}
\\addtolength{\\textwidth}{-2\\oddsidemargin}

% Uncommenting this gives more lines on each page but causes some
% "Underfull \\vbox" error messages.
%\\setlength{\\topmargin}{-1in}
%\\addtolength{\\topmargin}{2cm}
%\\setlength{\\textheight}{\\paperheight}
%\\addtolength{\\textheight}{-2in}
%\\addtolength{\\textheight}{-2\\topmargin}
%\\addtolength{\\textheight}{-\\headheight}
%\\addtolength{\\textheight}{-\\headsep}
%\\addtolength{\\textheight}{-\\footskip}
EOF

        if ( $opt_linenum ) {
            print $handle <<"EOF";

\\newcommand{\\linenum}[1]{\\makebox[0pt][r]{\\footnotesize #1}}
EOF
        }

        if ( $opt_pretty ) {
            my $time = localtime;
            print $handle <<"EOF";

\\lhead{$time}
\\chead{$filename}
\\rhead{Page~\\thepage}
\\lfoot{}
\\cfoot{}
\\rfoot{}

\\pagestyle{fancy}
EOF
        }
        print $handle <<"EOF";

\\begin{document}
EOF

    }
    print $handle <<"EOF";
{\\begin{tabbing}
EOF

}


########################################################################
# write_postamble
#
# Write the LaTeX document postamble.
#
sub write_postamble (*) {
    my $handle = shift;
    print $handle "\\end{tabbing}}\n";
    print $handle "\\end{document}\n" if $opt_standalone;
}


########################################################################
# split_chunk
#
# Split a chunk of Matlab code into a list where each element is a
# Matlab string ('...'), comment (%...) or neither.  Works on a single
# chunk of code or a list of code chunks.  If input is omitted, uses $_.
#
sub split_chunk (@) {

    my @ichunks  = @_ ? @_ : $_;        # Get input.
    my @ochunks = ();                   # Initialize output.

    foreach ( @ichunks ) {

        # Remove trailing whitespace (including newline character).
        s/\s+$//g;

        # Now do the splitting.
        my @parts                       # Each part may be
          = /   $ml_normal_regex        # neither comment nor string
              | $ml_string_regex        # or string
              | $ml_comment_regex       # or comment.
            /gox;

        # Append result to output list.
        push @ochunks, @parts;

    }

    # Return output argument(s).
    return          @ochunks if         wantarray;  # List context.
    return join "", @ochunks if defined wantarray;  # Scalar context.
    $_ =   join "", @ochunks;                       # Void context.

}

########################################################################
# process_line
#
# Converts a line of Matlab code to a line of LaTeX code.  If input is
# omitted, uses $_.
#
sub process_line (;$) {

    my $line = @_ ? $_[0] : $_;         # Get input.

    # Split chunk into Matlab strings, comments and the rest.
    my @parts = split_chunk $line;

    foreach my $part ( @parts ) {

        # The part is a string.
        if ( $part =~ /^'/ ) {                      #'
             $part = transl_tex_char $part;
             $part = $begin_string . $part . $end_string;
        }

        # This part is a comment.
        elsif ( $part =~ /^%/ ) {
            $part = transl_tex_char $part;
            $part = $begin_comment . $part . $end_comment;
        }

        # This part is neither string nor comment.
        else {

            #
            # Split the code into parts where each part contains at most
            # one keyword. The elements with keywords may also contain
            # non-word-characters surrouding the keyword.
            #
            my @sub_parts = split /($ml_keyword_regex)/gox, $part;

            my @new_parts = ();
            my $non_keyword = "";
            foreach my $sub_part ( @sub_parts ) {

                # Sometimes it is empty, but I don't know why...
                next if $sub_part eq "";

                # If this part does not contain a keyword, append it to
                # the chunk of code without keywords.
                if ( $sub_part !~ /$ml_keyword_regex/ox ) {
                    $non_keyword .= $sub_part;
                }

                # If it does contain a keyword, split this part into the
                # part before the keyword, the keyword, and the part
                # after the keyword.
                else {
                    my @pieces = split /\b/, $sub_part;
                    foreach my $piece ( @pieces ) {

                        # If this part is not a keyword, append it to
                        # the chunk of code without keywords.
                        if ( $piece !~ /^\w/ ) {
                            $non_keyword .= $piece;

                        # If this part is a keyword, convert the part
                        # before the keyword to TeX and prepend it to
                        # the list. Then prepend the keyword to the
                        # list.
                        } else {
                            push @new_parts,
                              $begin_normal,
                                transl_tex_char( $non_keyword ),
                                  $end_normal
                                    if $non_keyword ne "";
                            $non_keyword = "";
                            push @new_parts,
                              $begin_keyword,
                                $piece,
                                  $end_keyword;
                        }   # of if
                    }   # of foreach
                }   # of if
            }   # of foreach

            # There still might be some code after the last keyword that
            # we haven't processed, so do it now.
            push @new_parts, $begin_normal,
              transl_tex_char( $non_keyword ), $end_normal
                if $non_keyword ne "";

            # Now glue all the parts together.
            $part = join "", @new_parts;

        }   # of if
    }   # of foreach

    $line = join "", @parts;

    return $line if defined wantarray;
    $_ =   $line;

}

########################################################################
# process_code
#
# Convert Matlab code to LaTeX code.
#
sub process_code (**$) {

    my ( $ifh, $ofh, $filename ) = @_;

    write_preamble $ofh, $filename;     # Write preamble.

    if ( $opt_online ) {
        # Only the online help is of interest.
        while ( <$ifh> ) {
            last if /^%/;
        }
        while ( defined ) {
            last unless /^%/;
            process_line;
            print $ofh "\\linenum{$.}\\quad" if $opt_linenum;
            print $ofh "$_\\\\\n";
            $_ = <$ifh>;
        }
    } else {
        # The whole file should be converted.
        while ( <$ifh> ) {
            process_line;
            print $ofh "\\linenum{$.}\\quad" if $opt_linenum;
            print $ofh "$_\\\\\n";
        }

    }

    write_postamble $ofh;               # Write postamble.

}

########################################################################
# process_file
#
# Process a single Matlab m-file.
#
sub process_file ($) {

    # Get name of input file and try to open it.
    my $ifile = shift;
    my $ifh   = FileHandle->new($ifile)
      or warn( "$program: can't open $ifile: $!\n" ), return;

    #
    # Generate name of output file and try to open the file. Ignore case
    # on file suffix so it works properly on operating systems where
    # case doesn't matter.
    #
    ( my $ofile = $ifile ) =~ s/\.[mM]$/.tex/;
    if ( defined $outputdir ) {
        $ofile = basename $ofile;
        $ofile = "$outputdir$ofile";
    }

    my $ofh = FileHandle->new(">$ofile")
      or $ifh->close,
        warn( "$program: can't open $ofile: $!\n" ), return;

    # Print name of file we are about to process.
    print "$ifile -> $ofile\n" if $opt_verbose;

    process_code $ifh, $ofh, $ifile;

    $ifh->close;
    $ofh->close;

}

########################################################################
# process_files
#
# Processes all the Matlab m-files.
#
sub process_files () {

    if ( @ARGV ) {
        while ( my $file = shift @ARGV ) {
            # Perform file globbing if necessary.  This might be too simple.
            if ( $file =~ m! (^|/)~ | [*?{}[\]^] !x ) {
                my @files = glob $file;
                if ( @files ) {
                    $file = shift @files;       # Get first and prepend the
                    unshift @ARGV, @files;      # rest to the argument list.
                } else {
                    warn "$program: $file: No such file or directory\n";
                    next;
                }
            }
            process_file $file;
        }
    } else {
        process_code \*STDIN, \*STDOUT, "stdin";
    }

}

########################################################################
# define_regexes
#
# Define regular expressions used to match different parts of the Matlab
# code.
#
sub define_regexes () {

    #
    # Regex to match a chunk of Matlab source code which is neither a
    # Matlab string or comment. A single quote is a transpose operator
    # if it is preceded by a letter, digit, underscore, dot or closing
    # delimiter (parenthesis, brace or bracket).
    #
    $ml_normal_regex =
      q!
        (?:             # Non-backreferencing grouping parenthesis.
          [])}\w.]'+    # Either a character followed by one or more
                        #   single quotes that are transpose operators
            |           # or else
          [^'%]         #   any character except single quote (which
                        #   starts a string) or a percent sign (which
                        #   starts a comment).
        )+              # Match one or more times.
       !;

    #
    # Regex to match a Matlab string if it is known that the leading '
    # (single quote character) is not a transpose operator. Literal
    # single quote characters are made up by two consecutive single
    # quote characters, e.g., 'It''s great isn''t it?'.
    #
    $ml_string_regex =
      q!
        '               # Opening single quote that starts the string.
          [^'\n]*       # Zero or more chars that are neither single
                        #   quotes (special) nor newlines (illegal).
          (?:           # Non-backreferencing grouping parenthesis.
            ''          # An embedded (literal) single quote character.
            [^'\n]*     # Again, zero or more chars that are neither
                        #   single quotes nor newlines.
          )*            # Match zero or more times.
        '               # Closing single quote that ends the string.
       !;

    #
    # Regex to match a Matlab comment if it is known that the leading %
    # (percent sign) starts a comment. Comments continue to the end of
    # the line.
    #
    $ml_comment_regex =
      q!
        %               # Opening percent sign that starts the comment.
        [^\n]*          # Match to the end of the line.
       !;

    #
    # Regex to match a Matlab keyword and surrounding characters if it
    # is known that the target string does not contain Matlab strings.
    # When "end" is used as an array subscript it is not a keyword, but
    # it is not easy to get it right. The solution below catches some of
    # the cases, but not all.
    #
    $ml_keyword_regex =
      q!
          ^
          \s*
          function
        |
          (?:
              ^
            |
              [,;]
          )
          \s*
          (?:
              (?: break | return | switch | case
                | otherwise | if | else | elseif
                | for | while | try | catch
              )
              \b
            |
              end
              \s*
              (?:
                  [,;]
                |
                  $
              )
          )
       !;

}

########################################################################
# Main part.
########################################################################

read_rcfile;            # Read initialization file.
define_regexes;         # Define regular expressions.
process_files;          # Process the files.
