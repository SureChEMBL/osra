/*  GNU Ocrad - Optical Character Recognition program
    Copyright (C) 2003, 2004, 2005, 2006 Antonio Diaz Diaz.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*/
/*
    Return values: 0 for a normal exit, 1 for environmental problems
    (file not found, invalid flags, I/O errors, etc), 2 to indicate a
    corrupt or invalid input file, 3 for an internal consistency error
    (eg, bug) which caused ocrad to panic.
*/

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "arg_parser.h"
#include "common.h"
#include "rational.h"
#include "rectangle.h"
#include "page_image.h"
#include "textpage.h"


namespace {

// Date of this version: 2006-10-20

const char * invocation_name = 0;
const char * const Program_name    = "GNU Ocrad";
const char * const program_name    = "ocrad";
const char * const program_version = "0.16";
const char * const program_year    = "2006";

struct Input_control
  {
  Transformation transformation;
  int rindex, layout_level, scale;
  Rational threshold, ltrb[4];
  char specialtype;
  bool copy, crop, invert;

  Input_control() throw()
    : rindex( -1 ), layout_level( 0 ), scale( 0 ), threshold( -1 ),
      specialtype( 0 ), copy( false ), crop( false ), invert( false ) {}

  bool read_crop_rectangle( const char * ptr ) throw();
  bool set_threshold( const char * ptr ) throw();
  };


void show_error( const char * msg, const int errcode = 0, const bool help = false ) throw()
  {
  if( msg && msg[0] != 0 )
    {
    std::fprintf( stderr, "%s: %s", program_name, msg );
    if( errcode > 0 ) std::fprintf( stderr, ": %s", strerror( errcode ) );
    std::fprintf( stderr, "\n" );
    }
  if( help && invocation_name && invocation_name[0] != 0 )
    std::fprintf( stderr, "Try `%s --help' for more information.\n", invocation_name );
  }


bool Input_control::read_crop_rectangle( const char * ptr ) throw()
  {
  int c = ltrb[0].parse( ptr );				// left
  if( c && ltrb[0] >= 0 && ptr[c] == ',' )
    {
    int i = c + 1;
    c = ltrb[1].parse( &ptr[i] );			// top
    if( c && ltrb[1] >= 0 && ptr[i+c] == ',' )
      {
      i += c + 1; c = ltrb[2].parse( &ptr[i] );		// right
      if( c && ltrb[2] > 0 && ptr[i+c] == ',' )
        {
        i += c + 1; c = ltrb[3].parse( &ptr[i] );	// bottom
        if( c && ltrb[3] > 0 ) { crop = true; return true; }
        }
      }
    }
  show_error( "invalid crop rectangle", 0, true );
  return false;
  }


bool Input_control::set_threshold( const char * ptr ) throw()
  {
  Rational tmp;
  if( tmp.parse( ptr ) && tmp >= 0 && tmp <= 1 )
    { threshold = tmp; return true; }
  show_error( "threshold out of limits (0.0 - 1.0)", 0, true );
  return false;
  }


void show_help() throw()
  {
  std::printf( "%s - Optical Character Recognition program.\n", Program_name );
  std::printf( "Reads pnm file(s), or standard input, and sends text to standard output.\n" );
  std::printf( "\nUsage: %s [options] [files]\n", invocation_name );
  std::printf( "Options:\n" );
  std::printf( "  -h, --help               display this help and exit\n" );
  std::printf( "  -V, --version            output version information and exit\n" );
  std::printf( "  -a, --append             append text to output file\n" );
  std::printf( "  -b, --block=<n>          process only the specified text block\n" );
  std::printf( "  -c, --charset=<name>     try `--charset=help' for a list of names\n" );
  std::printf( "  -e, --filter=<name>      try `--filter=help' for a list of names\n" );
  std::printf( "  -f, --force              force overwrite of output file\n" );
  std::printf( "  -F, --format=<fmt>       output format (byte, utf8)\n" );
  std::printf( "  -i, --invert             invert image levels (white on black)\n" );
  std::printf( "  -l, --layout=<n>         layout analysis, 0=none, 1=column, 2=full\n" );
  std::printf( "  -o <file>                place the output into <file>\n" );
  std::printf( "  -p, --crop=<l,t,r,b>     crop input image by given rectangle\n" );
  std::printf( "  -s, --scale=[-]<n>       scale input image by [1/]<n>\n" );
  std::printf( "  -t, --transform=<name>   try `--transform=help' for a list of names\n" );
  std::printf( "  -T, --threshold=<n%%>     threshold for binarization (0-100%%)\n" );
  std::printf( "  -v, --verbose            be verbose\n" );
  std::printf( "  -x <file>                export OCR Results File to <file>\n" );
  if( Ocrad::verbose )
    {
    std::printf( "  -1..6                    pnm output file type (debug)\n" );
    std::printf( "  -C, --copy               'copy' input to output (debug)\n" );
    std::printf( "  -D, --debug=<level>      (0-100) output intermediate data (debug)\n" );
    std::printf( "  -S <type>                make a 'special file' (debug)\n" );
    }
  std::printf( "\nReport bugs to bug-ocrad@gnu.org\n" );
  }


void show_version() throw()
  {
  std::printf( "%s version %s\n", Program_name, program_version );
  std::printf( "Copyright (C) %s Antonio Diaz Diaz.\n", program_year );
  std::printf( "This program is free software; you may redistribute it under the terms of\n" );
  std::printf( "the GNU General Public License.  This program has absolutely no warranty.\n" );
  }


const char * my_basename( const char * filename ) throw()
  {
  const char * c = filename;
  while( *c ) { if( *c == '/' ) filename = c + 1; ++c; }
  return filename;
  }


int process_file( FILE *infile, const char * infile_name,
                  const Input_control & input_control,
                  const Control & control ) throw()
  {
  if( Ocrad::verbose )
    std::fprintf( stderr, "processing file `%s'\n", infile_name );
  try
    {
    Page_image page_image( infile, input_control.threshold, input_control.invert );

    if( input_control.crop &&
        !page_image.crop( input_control.ltrb, input_control.threshold ) )
      {
      if( Ocrad::verbose )
        std::fprintf( stderr, "file `%s' totally cropped out\n", infile_name );
      return 1;
      }
    page_image.transform( input_control.transformation );
    page_image.scale( input_control.scale );
    page_image.analyse_layout( input_control.layout_level );
    if( Ocrad::verbose )
      std::fprintf( stderr, "number of text blocks = %d\n", page_image.zones() );
    if( input_control.threshold < 0 ) page_image.adapt_thresholds();

    if( input_control.rindex >= page_image.zones() )
      {
      std::fprintf( stderr,"This page has only %d text block(s)\n", page_image.zones() );
      return 1;
      }
    if( input_control.specialtype != 0 )
      {
      if( control.outfile )
        {
        const char & t = input_control.specialtype;
        if( t == 'v' || t == 'h' ) page_image.histogramize( t == 'v' );
        else { show_error( "bad special type" ); return 1; }
        page_image.save( control.outfile, control.filetype );
        }
      return 0;
      }
    if( input_control.copy )
      {
      if( control.outfile )
        {
        if( input_control.layout_level == 0 && page_image.zones() == 1 )
          page_image.save( control.outfile, control.filetype );
        else
          for( int c = 0; c < page_image.zones(); ++c )
            if( input_control.rindex < 0 || input_control.rindex == c )
              page_image.save( control.outfile, control.filetype, c );
        }
      return 0;
      }

    Textpage textpage( page_image, my_basename( infile_name ), control );
    }
  catch( Page_image::Error e ) { show_error( e.s ); return 2; }
  if( Ocrad::verbose ) std::fprintf( stderr, "\n" );
  return 0;
  }

} // end namespace


bool Ocrad::verbose = false;

void Ocrad::internal_error( const char * msg ) throw()
  {
  std::string s( "internal error: " ); s += msg;
  show_error( s.c_str() );
  exit( 3 );
  }


// 'infile' contains the scanned image (in pnm format) to be converted
// to text.
// 'outfile' is the destination for the text version of the scanned
// image. (or for a pnm file if debugging).
// 'exportfile' is the Ocr Results File.
//
/*int main( int argc, char * argv[] ) throw()
  {
  Input_control input_control;
  Control control;
  const char *outfile_name = 0, *exportfile_name = 0;
  bool append = false, force = false;
  invocation_name = argv[0];

  static const Arg_parser::Option options[] =
    {
    { '1', 0,           Arg_parser::no  },
    { '2', 0,           Arg_parser::no  },
    { '3', 0,           Arg_parser::no  },
    { '4', 0,           Arg_parser::no  },
    { '5', 0,           Arg_parser::no  },
    { '6', 0,           Arg_parser::no  },
    { 'a', "append",    Arg_parser::no  },
    { 'b', "block",     Arg_parser::yes },
    { 'c', "charset",   Arg_parser::yes },
    { 'C', "copy",      Arg_parser::no  },
    { 'D', "debug",     Arg_parser::yes },
    { 'e', "filter",    Arg_parser::yes },
    { 'f', "force",     Arg_parser::no  },
    { 'F', "format",    Arg_parser::yes },
    { 'h', "help",      Arg_parser::no  },
    { 'i', "invert",    Arg_parser::no  },
    { 'l', "layout",    Arg_parser::yes },
    { 'o', 0,           Arg_parser::yes },
    { 'p', "crop",      Arg_parser::yes },
    { 's', "scale",     Arg_parser::yes },
    { 'S', 0,           Arg_parser::yes },
    { 't', "transform", Arg_parser::yes },
    { 'T', "threshold", Arg_parser::yes },
    { 'v', "verbose",   Arg_parser::no  },
    { 'V', "version",   Arg_parser::no  },
    { 'x', 0,           Arg_parser::yes },
    {  0 , 0,           Arg_parser::no  } };

  Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  int argind;
  for( argind = 0; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;					// no more options
    const char * arg = parser.argument( argind ).c_str();
    switch( code )
      {
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6': control.filetype = code; break;
      case 'C': input_control.copy = true; break;
      case 'D': control.debug_level = std::strtol( arg, 0, 0 ); break;
      case 'F': if( !control.set_format( arg ) )
                  { show_error( "bad output format", 0, true ); return 1; }
                break;
      case 'S': input_control.specialtype = arg[0]; break;
      case 'T': if( !input_control.set_threshold( arg ) ) return 1; break;
      case 'V':	show_version(); return 0;
      case 'a': append = true; break;
      case 'b': input_control.rindex = std::strtol( arg, 0, 0 ) - 1; break;
      case 'c': if( !control.charset.enable( arg ) )
                  { control.charset.show_error( program_name, arg ); return 1; }
                break;
      case 'e': if( !control.filter.set( arg ) )
                  { control.filter.show_error( program_name, arg ); return 1; }
                break;
      case 'f': force = true; break;
      case 'h': show_help(); return 0;
      case 'i': input_control.invert = true; break;
      case 'l': input_control.layout_level = std::strtol( arg, 0, 0 ); break;
      case 'o':	outfile_name = arg; break;
      case 'p': if( !input_control.read_crop_rectangle( arg ) ) return 1; break;
      case 's': input_control.scale = std::strtol( arg, 0, 0 ); break;
      case 't': if( !input_control.transformation.set( arg ) )
                  { input_control.transformation.show_error( program_name, arg );
                  return 1; }
                break;
      case 'v': Ocrad::verbose = true; break;
      case 'x':	exportfile_name = arg; break;
      default : Ocrad::internal_error( "uncaught option" );
      }
    } // end process options

  if( outfile_name && std::strcmp( outfile_name, "-" ) != 0 )
    {
    if( append ) control.outfile = std::fopen( outfile_name, "a" );
    else if( force ) control.outfile = std::fopen( outfile_name, "w" );
    else if( ( control.outfile = std::fopen( outfile_name, "wx" ) ) == 0 )
      {
      std::fprintf( stderr, "Output file %s already exists.\n", outfile_name );
      return 1;
      }
    if( !control.outfile )
      { std::fprintf( stderr, "Cannot open %s\n", outfile_name ); return 1; }
    }

  if( exportfile_name && control.debug_level == 0 &&
      input_control.scale == 0 && input_control.specialtype == 0 &&
      !input_control.copy )
    {
    if( std::strcmp( exportfile_name, "-" ) == 0 )
      { control.exportfile = stdout; if( !outfile_name ) control.outfile = 0; }
    else
      {
      control.exportfile = std::fopen( exportfile_name, "w" );
      if( !control.exportfile )
        {
        std::fprintf( stderr, "Cannot open %s\n", exportfile_name );
        return 1;
        }
      }
    std::fprintf( control.exportfile,
                  "# Ocr Results File. Created by %s version %s\n",
                  Program_name, program_version );
    }

  // process any remaining command line arguments (input files)
  FILE *infile = (argind < parser.arguments()) ? 0 : stdin;
  const char *infile_name = "-";
  int retval = 0;
  while( true )
    {
    while( infile != stdin )
      {
      if( infile ) std::fclose( infile );
      if( argind >= parser.arguments() ) { infile = 0; break; }
      infile_name = parser.argument( argind++ ).c_str();
      if( std::strcmp( infile_name, "-" ) == 0 ) infile = stdin;
      else infile = std::fopen( infile_name, "r" );
      if( infile ) break;
      std::fprintf( stderr, "Cannot open %s\n", infile_name );
      if( retval == 0 ) retval = 1;
      }
    if( !infile ) break;

    int tmp = process_file( infile, infile_name, input_control, control );
    if( infile == stdin )
      {
      if( tmp <= 1 )
        {
        int ch;
        do ch = std::fgetc( infile ); while( std::isspace( ch ) );
        std::ungetc( ch, infile );
        }
      if( tmp > 1 || std::feof( infile ) || std::ferror( infile ) ) infile = 0;
      }
    if( tmp > retval ) retval = tmp;
    if( control.outfile ) std::fflush( control.outfile );
    if( control.exportfile ) std::fflush( control.exportfile );
    }
  if( control.outfile ) std::fclose( control.outfile );
  if( control.exportfile ) std::fclose( control.exportfile );
  return retval;
  }
*/
