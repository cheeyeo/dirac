/* ***** BEGIN LICENSE BLOCK *****
*
* $Id: YUV420Down2x2.cpp,v 1.3 2008/10/01 00:50:27 asuraparaju Exp $ $Name: Dirac_1_0_2 $
*
* Version: MPL 1.1/GPL 2.0/LGPL 2.1
*
* The contents of this file are subject to the Mozilla Public License
* Version 1.1 (the "License"); you may not use this file except in compliance
* with the License. You may obtain a copy of the License at
* http://www.mozilla.org/MPL/
*
* Software distributed under the License is distributed on an "AS IS" basis,
* WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
* the specific language governing rights and limitations under the License.
*
* The Original Code is BBC Research and Development code.
*
* The Initial Developer of the Original Code is the British Broadcasting
* Corporation.
* Portions created by the Initial Developer are Copyright (C) 2008.
* All Rights Reserved.
*
* Contributor(s):
*
* Alternatively, the contents of this file may be used under the terms of
* the GNU General Public License Version 2 (the "GPL"), or the GNU Lesser
* Public License Version 2.1 (the "LGPL"), in which case the provisions of
* the GPL or the LGPL are applicable instead of those above. If you wish to
* allow use of your version of this file only under the terms of the either
* the GPL or LGPL and not to allow others to use your version of this file
* under the MPL, indicate your decision by deleting the provisions above
* and replace them with the notice and other provisions required by the GPL
* or LGPL. If you do not delete the provisions above, a recipient may use
* your version of this file under the terms of any one of the MPL, the GPL
* or the LGPL.
* ***** END LICENSE BLOCK ***** */

/*****************************************************************
File YUV420Down2x2.cpp

Utility for filtering a sequence of frames, stored in a single
file in raw YUV420 format, to a single output file with half dimensions,
a (-1,9,9,-1) filter and subsampling.
This utility is a filter taking input on stdin and generating its
output on stdout.

Original author: Thomas Davies (based on filters by Tim Borer)
****************************************************************/
#include <libdirac_common/arrays.h>
#include <cmath>
#include <stdlib.h> //Contains EXIT_SUCCESS, EXIT_FAILURE
#include <iostream> //For cin, cout, cerr
#include <algorithm> //For fill_n
using std::cout;
using std::cin;
using std::cerr;
using std::clog;
using std::endl;
using std::ios_base;
using std::fill_n;

#include <setstdiomode.h>
using namespace dirac_vu;
using namespace dirac;

void VFilter( const TwoDArray<unsigned char>& pic_data, TwoDArray<unsigned char>& out_data, const OneDArray<int>& filter, const int bits );
void HFilter( const TwoDArray<unsigned char>& pic_data, TwoDArray<unsigned char>& out_data, const OneDArray<int>& filter, const int bits );

double sinxoverx( const double val )
{
    if ( 0.0f == val )
        return 1.0;
    else
        return sin(val)/val;

}

OneDArray<int> MakeLPRectFilter( const float bw, const int bits )
{
    const int tl = 8;
    const float pi = 3.1415926535;

    OneDArray<double> double_filter( Range( -tl, tl ) );
    OneDArray<int> int_filter( Range( -tl, tl) );

    // Use the Hanning window
    for (int i=double_filter.First(); i<=double_filter.Last(); ++i)
    {
        double_filter[i] = cos( (pi*i)/
                                   (double_filter.Length()+1) );
    }

    // Apply sinc function
    for (int i=double_filter.First(); i<=double_filter.Last(); ++i)
    {
        double_filter[i] *= sinxoverx( pi*1.0*bw*i );
    }

    // Get DC gain = 1<<bits
    double sum = 0.0;
    for (int i=double_filter.First(); i<=double_filter.Last(); ++i)
        sum += double_filter[i];

    for (int i=double_filter.First(); i<=double_filter.Last(); ++i)
    {
        double_filter[i] *= double(1<<(bits+4));
        double_filter[i] /= sum;
    }

    // Turn the float filter into an integer filter
    for (int i=double_filter.First(); i<=double_filter.Last(); ++i)
    {
        int_filter[i] = double_filter[i]>0 ? int( double_filter[i]+0.5 ) : -int( -double_filter[i]+0.5 );
        int_filter[i] = (int_filter[i]+8)>>4;
    }

    return int_filter;
}

void HFilter( const TwoDArray<unsigned char>& pic_data, TwoDArray<unsigned char>& out_data, const OneDArray<int>& filter, const int bits )
{
    unsigned char* line_data;
    const int offset = (1<<(bits-1));

    int sum;
    int i, iout;

    for (int j=0; j<pic_data.LengthY(); ++j)
    {
        line_data = out_data[j];
        // Do the first bit
        for (i=0, iout=0; i<filter.Last(); i+=2, iout++)
        {
            sum = offset;
            for (int k=filter.Last(); k>=filter.First(); --k)
                sum += filter[k]*static_cast<int>(pic_data[j][std::max(i-k,0)]);
            sum >>= bits;
            sum = std::min( 255, std::max( 0, sum) );
            line_data[iout] = static_cast<unsigned char>( sum );
        }// i
        // Do the middle bit
        for (; i<=pic_data.LastX()+filter.First(); i+=2, iout++)
        {
            sum = offset;
            for (int k=filter.Last(); k>=filter.First(); --k)
                sum += filter[k]*static_cast<int>( pic_data[j][i-k] );
            sum >>= bits;
            sum = std::min( 255, std::max( 0, sum) );
            line_data[iout] = static_cast<unsigned char>( sum );
        }// i
        // Do the last bit
	for (; i<pic_data.LengthX(); i+=2, iout++)
	{
            sum = offset;
            for (int k=filter.Last(); k>=filter.First(); --k)
                sum += filter[k]*static_cast<int>( pic_data[j][std::min(i-k,pic_data.LastX())] );
            sum >>= bits;
            sum = std::min( 255, std::max( 0, sum) );
            line_data[iout] = static_cast<unsigned char>( sum );
        }// i
    }// j

}

void VFilter( const TwoDArray<unsigned char>& pic_data, TwoDArray<unsigned char>& out_data, const OneDArray<int>& filter, const int bits )
{
    const int offset = (1<<(bits-1));

    int sum, j, jout;

    // Do the first bit
    for (j=0, jout=0; j<filter.Last(); j+=2, jout++)
    {

        for (int i=0; i<pic_data.LengthX(); ++i)
        {
            sum = offset;
            for (int k=filter.Last(); k>=filter.First(); --k)
                sum += filter[k]*static_cast<int>( pic_data[std::max(j-k,0)][i] );
            sum >>= bits;
            sum = std::min( 255, std::max( 0, sum) );
            out_data[jout][i] = static_cast<unsigned char>( sum );
        }// i

    }// j

    // Do the middle bit
    for (; j<=pic_data.LastY()+filter.First(); j+=2, jout++)
    {

        for (int i=0; i<pic_data.LengthX(); ++i)
        {
            sum = offset;
            for (int k=filter.Last(); k>=filter.First(); --k)
                sum += filter[k]*static_cast<int>( pic_data[j-k][i] );
            sum >>= bits;
            sum = std::min( 255, std::max( 0, sum) );
            out_data[jout][i] = static_cast<unsigned char>( sum );
        }// i

    }// j

    // Do the last bit
    for (; j<pic_data.LengthY(); j+=2, jout++)
    {

        for (int i=0; i<pic_data.LengthX(); ++i)
        {
            sum = offset;
            for (int k=filter.Last(); k>=filter.First(); --k)
                sum += filter[k]*static_cast<int>( pic_data[std::min(j-k,pic_data.LastY())][i] );
            sum >>= bits;
            sum = std::min( 255, std::max( 0, sum) );
            out_data[jout][i] = static_cast<unsigned char>( sum );
        }// i

    }// j
}


void h_filter(unsigned char *in_array, unsigned char *out_array, const int w, const int h);

void v_filter(unsigned char *in_array, unsigned char *out_array, const int w, const int h);

int main(int argc, char * argv[] ) {

    if (argc != 4) {
        cout << "\"YUV420Down2x2\" command line format is:" << endl;
        cout << "    Argument 1: width (pixels) e.g. 720" << endl;
        cout << "    Argument 2: height (lines) e.g. 576" << endl;
        cout << "    Argument 3: number of frames e.g. 3" << endl;
        cout << "    Example: YUV420Down2x2 <foo >bar 720 576 3" << endl;
        cout << "        converts 3 frames, of 720x576 pixels, from file foo to file bar" << endl;
        return EXIT_SUCCESS; }

    //Get command line arguments
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);
    int frames = atoi(argv[3]);

    //Set standard input and standard output to binary mode.
    //Only relevant for Windows (*nix is always binary)
    if ( setstdinmode(std::ios_base::binary) == -1 ) {
        cerr << "Error: could not set standard input to binary mode" << endl;
        return EXIT_FAILURE; }
    if ( setstdoutmode(std::ios_base::binary) == -1 ) {
        cerr << "Error: could not set standard output to binary mode" << endl;
        return EXIT_FAILURE; }

    //Allocate memory for input and output buffers.
    const int ybuffersizein = height*width;
    TwoDArray<unsigned char> ybufferin( height, width );
    const int uvbuffersizein = height*width/4;
    TwoDArray<unsigned char> ubufferin ( height/2, width/2);
    TwoDArray<unsigned char> vbufferin ( height/2, width/2);

    int width2 = width/2;
    int height2 = height/2;

    // Output buffers
    const int ybuffersizeout = height2*width2;
    TwoDArray<unsigned char> ybufferout( height2, width2);
    const int uvbuffersizeout = height2*width2/4;
    TwoDArray<unsigned char> ubufferout(height2/2, width2/2);
    TwoDArray<unsigned char> vbufferout(height2/2, width2/2);

    // Temporary buffers
    TwoDArray<unsigned char> ybuffertemp( height, width2);
    TwoDArray<unsigned char> uvbuffertemp(height/2, width2/2);

    // filter with 14-bit accuracy
    OneDArray<int> filter=MakeLPRectFilter(0.5, 16);

    //Create references for input and output stream buffers.
    //IO is via stream buffers for efficiency
    std::streambuf& inbuf = *(cin.rdbuf());
    std::streambuf& outbuf = *(cout.rdbuf());

    for (int frame=0; frame<frames; ++frame) {

        clog << "Processing frame " << (frame+1) << "\r";

        //Read frames of Y then U then V components
        if ( inbuf.sgetn(reinterpret_cast<char*>(&ybufferin[0][0]), ybuffersizein) < ybuffersizein ) {
            cerr << "Error: failed to read Y component of frame " << frame << endl;
            return EXIT_FAILURE; }
        if ( inbuf.sgetn(reinterpret_cast<char*>(&ubufferin[0][0]), uvbuffersizein) < uvbuffersizein ) {
            cerr << "Error: failed to read U component of frame " << frame << endl;
            return EXIT_FAILURE; }
        if ( inbuf.sgetn(reinterpret_cast<char*>(&vbufferin[0][0]), uvbuffersizein) < uvbuffersizein ) {
            cerr << "Error: failed to read V component of frame " << frame << endl;
            return EXIT_FAILURE; }

        HFilter( ybufferin, ybuffertemp, filter, 16 );
        VFilter( ybuffertemp, ybufferout, filter, 16 );

        HFilter( ubufferin, uvbuffertemp, filter, 16 );
        VFilter( uvbuffertemp, ubufferout, filter, 16 );

        HFilter( vbufferin, uvbuffertemp, filter, 16 );
        VFilter( uvbuffertemp, vbufferout, filter, 16 );

        //Write frames of YUV
        if ( outbuf.sputn(reinterpret_cast<char*>(&ybufferout[0][0]), ybuffersizeout) < ybuffersizeout ) {
            cerr << "Error: failed to write frame " << frame << endl;
            return EXIT_FAILURE; }
        if ( outbuf.sputn(reinterpret_cast<char*>(&ubufferout[0][0]), uvbuffersizeout) < uvbuffersizeout ) {
            cerr << "Error: failed to write frame " << frame << endl;
            return EXIT_FAILURE; }
        if ( outbuf.sputn(reinterpret_cast<char*>(&vbufferout[0][0]), uvbuffersizeout) < uvbuffersizeout ) {
            cerr << "Error: failed to write frame " << frame << endl;
            return EXIT_FAILURE; }
    }

    return EXIT_SUCCESS;
}

void v_filter(unsigned char *in_array, unsigned char *out_array, const int w, const int h){

    int height2 = h/2;

    int val;

    // top line
    for (int x=0; x<w; ++x ){
        val = (-in_array[x]+9*in_array[x]+9*in_array[w+x]-
                         in_array[2*w+x]+8)>>4;
	val = std::min( 255, std::max(0,val ) );
	out_array[x] = static_cast<unsigned char>( val );

    }
    // middle lines
    for (int line_out=1,line_in=2; line_out<height2-1; line_out++,line_in+=2) {
        for (int x=0; x<w; ++x ){
            val = (-in_array[(line_in-1)*w+x]+
                   9*in_array[line_in*w+x]+
                   9*in_array[(line_in+1)*w+x]-
                   in_array[(line_in+2)*w+x]+8)>>4;
	    val = std::min( 255, std::max(0,val ) );
	    out_array[line_out*w+x] = static_cast<unsigned char>( val );
        }
    }

    // bottom line
    for (int x=0; x<w; ++x ){
        val = (-in_array[(h-3)*w+x]+
               9*in_array[(h-2)*w+x]+
               9*in_array[(h-1)*w+x]-
               in_array[(h-1)*w+x]+8)>>4;
	val = std::min( 255, std::max(0,val ) );
	out_array[(height2-1)*w+x] = static_cast<unsigned char>( val );
    }
}

void h_filter(unsigned char *in_array, unsigned char *out_array, const int w, const int h){

    int width2 = w/2;
    int val;

    for (int j=0; j<h; ++j) {

        val = (-in_array[j*w]+9*in_array[j*w]+9*in_array[j*w+1]-
                         in_array[j*w+2]+8)>>4;
	val = std::min( 255, std::max(0,val ) );
	out_array[j*width2] = static_cast<unsigned char>( val );

        for (int xpos_out=1, xpos_in=2; xpos_out<width2-1; xpos_out++, xpos_in+=2 ){
            val = (-in_array[j*w+xpos_in-1]+
                   9*in_array[j*w+xpos_in]+
                   9*in_array[j*w+xpos_in+1]-
                   in_array[j*w+xpos_in+2]+8)>>4;
	    val = std::min( 255, std::max(0,val ) );
	    out_array[j*width2+xpos_out] = static_cast<unsigned char>( val );
        }

        val = (-in_array[j*w+w-3]+
                9*in_array[j*w+w-2]+
                9*in_array[j*w+w-1]-
                in_array[j*w+w-1]+8)>>4;
        val = std::min( 255, std::max(0,val ) );
	out_array[j*width2+width2-1] = static_cast<unsigned char>( val );
    }
}
