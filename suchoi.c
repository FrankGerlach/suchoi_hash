/* Copyright (c) 2016, Frank Gerlach ( Frank_Gerlach@epam.com )
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors.
*/


/*******************************************************
 * The High Performance Hash Function "suchoi".
 * 
 * This function will hash ANY data structure so that
 * 
 * A) hash function runtime is acceptable
 * 
 * B) the hash output has very good distribution
 *    NOTE: Linear Congruential Generators do NOT have
 *          this property and will often yield highly 
 *          unbalanced Hash Tables !
 *
 * C) Both runtime and output quality is better than
 *    Adler32, especially on tiny processors 
 *    without hardware divison/modulo logic.
 * 
 * 
 * Note for users: On "large" computers, the s-box
 *      could be made bigger (256 elements) and also
 *      64 bit state could be used.
 *      I do not yet have experimental data how this
 *      affects performance.
 * 
 * 
 * Author: Frank Gerlach, EPAM LLPD, Minsk
 *
 ******************************************************/


#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

typedef uint32_t (*hashFuncType)(const char*,uint32_t);


// the first digits of PI, use as an s-box (see DES Standard to see what an s-box is)
//  (in short: a highly non-linear function)
const uint32_t c_pi[]=
{
    0x243F6A88,
    0x85A308D3,
    0x3198A2E0,
    0x3707344A,
    0x40938222,
    0x99F31D00,
    0x82EFA98E,
    0xC4E6C894,
    0x52821E63,
    0x8D01377B,
    0xE5466CF3,
    0x4E90C6CC,
    0x0AC29B7C,
    0x97C50DD3,
    0xF84D5B5B,
    0x54709179 };

uint32_t shapashnik(uint32_t input,uint32_t key);
    
/* efficient hash function. See head of file for details 
 * 
 * Explanation of "inner workings":
 * The S-box, as wide as the state, ensures that a single
 * bit change in the input will on average flip half the
 * bits of the output ("total avalanche effect").
 * 
 * Rotating the state will make sure that successive
 * identical input octets will not cancel each other out.
 * 
 * After 16 octets of input, rotation of the state can no
 * longer ensure that identical input octets will not cancel 
 * each other out. Therefore, we apply the shapashnik permutation
 * function on the state.
 * 
 * (For a 64 bit state, self-hashing would be needed after
 *  32 octets of input)
 *
 * Named after Pavel Suchoi, engineer.
 * */
uint32_t suchoi(const char* input, uint32_t input_size)
{
    
    uint32_t state = 0;
    uint32_t i;
    for(i=0; i < input_size; i++)
    {
        char octet = input[i];// ^ (uint8_t)state;
        uint8_t upperNibble = octet >>4;
        uint8_t lowerNibble = octet & 0xF;
        state ^= c_pi[upperNibble];
        state = (state << 31) | (state >> 1);//rotate state
        state ^= c_pi[lowerNibble];
        state = (state << 31) | (state >> 1);//rotate state
        if( (i & 0xf) == 0xf )//danger of "xor-cancellation" -> permute the state
        {
            state = shapashnik(state,state);             
        }
    }
    return state;
}

void int2octets(uint32_t input,uint8_t* output)
{
   uint8_t i=0;
   for(i=0; i < 4;i++)
   {
       output[i] = input & 0xFF;
       input >>= 8;
   }
}

uint32_t octets2int(uint8_t* input)
{
   uint32_t output = 0;
   int8_t i;
   for(i=3; i > 0 ;i--)
   {
       output |= input[i];
       output <<= 8;
   }
   output |= input[i];
   return output;
}

void printBin(uint32_t input)
{
   int32_t i;
   for(i=31; i >=0; i--)
   {
      if( input & (1 << i) )
      {
        printf("1");
      }
      else
      {
        printf("0");
      }
      if( i && ((i & 0x7) == 0) )
      {
         printf(".");
      }
   }
   printf("\n");
}


uint32_t rotate(uint32_t input,uint8_t count)
{
   return (input << count) | (input >> (32 - count));
}

uint32_t  swapBits(uint32_t input, uint8_t pos1, uint8_t pos2)
{
   uint32_t stencil1 = ( 1 << pos1 );
   uint32_t stencil2 = ( 1 << pos2 );
   uint32_t v1 = input & stencil1;
   uint32_t v2 = input & stencil2;
   uint32_t output = input;
   output &= (~stencil1);
   output &= (~stencil2);
   if( v2 )
   {
     output |= stencil1;
   }
   if( v1 )
   {
     output |= stencil2;
   }
   return output;
}

uint32_t shapashnik_lp(uint32_t input,uint32_t key)
{
  uint32_t output = 0;
  uint8_t i;
  for(i=0; i < 32;i++)
  {
    uint8_t index = key & 0x1F;
    key = rotate(key,5);
  
    swapBits(input,i,index); 
  }
  return output;      
} 

/*high performance permutation/diffusion function 

 Named after Barys Shapashnik, engineer. 
*/
uint32_t shapashnik(uint32_t input,uint32_t key)
{
   //printf("shapashnik i:%x\n",input);
   uint8_t i;
   uint32_t output;
   for(i=0; i < 10;i++)
   {
     uint32_t upper = input >> 16;
     uint32_t lower = input & 0xFFFF;
     output = (lower << 16) | upper;
     output = rotate(output,key & 0x7);
     key = rotate(key,3); 
     input = output;
   }
   //printf("shapashnik o:%x\n",output);
   return output;
}


/* calculate the number of different bits in a and b */
uint8_t bit_diff(uint32_t a,uint32_t b)
{
    uint8_t i; 
    uint8_t diffCtr=0;
    for(i=0;i < 32;i++)
    {
        if( (a & 1) == (b & 1) )
        {
           diffCtr++; 
        }
        a >>= 1;
        b >>= 1;
    }
    return diffCtr;
}





#ifdef UNITTEST_SUPERHASH

#include "../esp8266-wakaama/string/llpd_string.h"
#include <stdlib.h>
#include <sys/time.h>

/* the Adler32 function as a benchmark competitor */
const int MOD_ADLER = 65521;

uint32_t adler32(const char *data, uint32_t len) /* where data is the location of the data in physical memory and 
                                                      len is the length of the data in bytes */
{
    uint32_t a = 1, b = 0;
    size_t index;
    
    /* Process each byte of the data in order */
    for (index = 0; index < len; ++index)
    {
        a = (a + data[index]) % MOD_ADLER;
        b = (b + a) % MOD_ADLER;
    }
    
    return (b << 16) | a;
}

uint64_t getTimeStamp()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t ret = tv.tv_sec;
    ret *= 1000000;
    ret += tv.tv_usec;
    return ret;
}

void testHashFuncForHashtable( hashFuncType hashFunc )
{
    //perform an analysis of performance for use in a hashtable  
    ALLOC_STR(inputL3,0)
    
    gbAddStr(&inputL3,"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"); 
    
    const int numTests = 200000;
    const int tableSize = numTests * 2;
    uint8_t* ctrField = malloc(tableSize);
    if( ctrField == NULL)
    {
        return;
    }
    uint32_t i;
    for(i=0; i < tableSize; i++)
    {
        ctrField[i] =0;
    }
    
    /*
    for(i=0; i < (numTests/2); i++)
    {
        inputL3.rawBuf[16] = i & 0xFF;
        inputL3.rawBuf[16+1] = i >> 8;
        uint32_t hv = hashFunc(inputL3.rawBuf,inputL3.valid);
        uint32_t index = hv % tableSize;
        if(ctrField[index] != 255)
        {
            ctrField[index]++;
        }
    }
    
    for(i=0; i < (numTests/4); i++)
    {
        inputL3.rawBuf[32] = i & 0xFF;
        inputL3.rawBuf[55+1] = i >> 8;
        uint32_t hv = hashFunc(inputL3.rawBuf,inputL3.valid);
        uint32_t index = hv % tableSize;
        if(ctrField[index] != 255)
        {
            ctrField[index]++;
        }
    }*/
    
    getTimeStamp();
    
    uint64_t start = getTimeStamp();
    //printf("Seconds since Jan. 1, 1970: %ld\n", tv.tv_sec);
    
    //random string, random length
    srand(76273);
    
    for(i=0; i < (numTests); i++)
    {
        clearGoodBuf(&inputL3);
        uint32_t j=0;
        uint32_t count = rand();
        if( count < 0)
        {
            count *= -1;
        }
        count %= 300;
        for(j=0; j < count;j++)
        {
            gbAddChar(&inputL3,(char)rand());
        }
        
        
        uint32_t hv = hashFunc(inputL3.rawBuf,inputL3.valid);
        uint32_t index = hv % tableSize;
        if(ctrField[index] != 255)
        {
            ctrField[index]++;
        }
    }
    
    uint32_t histogramm[256];
    for(i=0; i < 256; i++)
    {
        histogramm[i]=0;
    }
    for(i=0; i < tableSize; i++)
    {
        histogramm[ctrField[i]]++;
    }
    
    printf("histogramm of bin list lengths:\n");
    for(i=0; i < 20/*256*/; i++)
    {
        printf("list length %d: %d\n",i,histogramm[i]);
    }
    
    freeGoodBuf(&inputL3);
    
    uint64_t end = getTimeStamp();
    
    printf("runtime: %d us\n", ((uint32_t)(end-start)));
    
    free(ctrField);
}

int main(int argc,char** argv)
{
    uint8_t buffer[4];
    int2octets(0xABCD0123,buffer); 
    uint8_t i=0;
    assert( (buffer[i++]== 0x23) );
    assert( (buffer[i++]== 0x01) );
    assert( (buffer[i++]== 0xcd) );
    assert( (buffer[i++]== 0xAB) );

    assert( octets2int(buffer) == 0xABCD0123);


    shapashnik(0x12345678,0x12345678);

    const char* input = "hallo welt";
    printf("t1:%x\n",suchoi(input,strlen(input)));
    
    const char* input2 = "hbllo welt";
    printf("t2:%x\n",suchoi(input2,strlen(input2)));
    
    const char* inputL1 = "Historically, poor choices had led to ineffective implementations of LCGs. A particularly illustrative example of this is RANDU, which was widely used in the early 1970s and led to many results which are currently being questioned because of the use of this poor LCG.[4]";
    printf("t2-1:%x\n",suchoi(inputL1,strlen(inputL1)));
    
    const char* inputL2 = "Historicallz, poor choices had led to ineffective implementations of LCGs. A particularly illustrative example of this is RANDU, which was widely used in the early 1970s and led to many results which are currently being questioned because of the use of this poor LCG.[4]";
    printf("t3-2:%x\n",suchoi(inputL2,strlen(inputL2)));

    ALLOC_STR(inputL3,0)
    gbAddStr(&inputL3,"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"); 
    uint32_t hv1 = suchoi(inputL3.rawBuf,inputL3.valid);
    
    
    
    
    inputL3.rawBuf[0] ^= 1;//change one bit
    uint32_t hv2 = suchoi(inputL3.rawBuf,inputL3.valid);
    uint8_t numDiff = bit_diff(hv1,hv2);
    printf("t4: diff: %d\n",numDiff);
    
    inputL3.rawBuf[1] ^= 1;//change one bit
    uint32_t hv3 = suchoi(inputL3.rawBuf,inputL3.valid);
    
    printf("t4: diff 2: %d\n",bit_diff(hv1,hv3));
    printf("t4: diff 3: %d\n",bit_diff(hv2,hv3));
    
    freeGoodBuf(&inputL3);
    
    printf("suchoi:\n");
    testHashFuncForHashtable(suchoi);
    printf("adler32:\n");
    testHashFuncForHashtable(adler32);

    
         
    
    return 1;
}


#endif
