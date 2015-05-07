/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be
 *  found  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

////////////////////////////////////////////////////////////////////
/*
 *  The code for FAST in this file is taken from
 *  http://www.edwardrosten.com/work/fast.html with the following
 *  copyright notice:
 *
 *  Copyright (c) 2006, 2008, 2009, 2010 Edward Rosten
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *  *Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  *Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  *Neither the name of the University of Cambridge nor the names of
 *  its contributors may be used to endorse or promote products derived
 *  from this software without specific prior written permission.

 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "vp9_corner_detect.h"

typedef struct {
  int x, y;
} xy;
typedef unsigned char byte;

xy* fast_nonmax(const byte* im, int xsize, int ysize, int stride,
                xy* corners, int numcorners, int barrier, int* numnx);
xy* fast_corner_detect_9(const byte* im, int xsize, int ysize, int stride,
                         int barrier, int* numcorners);

xy* fast_corner_detect_9(const byte* im, int xsize, int ysize, int stride,
                         int barrier, int* num) {
  int boundary = 3, y, cb, c_b;
  const byte  *line_max, *line_min;
  int rsize=512, total=0;
  xy *ret = (xy*)malloc(rsize*sizeof(xy));
  const byte* cache_0;
  const byte* cache_1;
  const byte* cache_2;
  int pixel[16];
  pixel[0] = 0 + 3 * stride;
  pixel[1] = 1 + 3 * stride;
  pixel[2] = 2 + 2 * stride;
  pixel[3] = 3 + 1 * stride;
  pixel[4] = 3 + 0 * stride;
  pixel[5] = 3 + -1 * stride;
  pixel[6] = 2 + -2 * stride;
  pixel[7] = 1 + -3 * stride;
  pixel[8] = 0 + -3 * stride;
  pixel[9] = -1 + -3 * stride;
  pixel[10] = -2 + -2 * stride;
  pixel[11] = -3 + -1 * stride;
  pixel[12] = -3 + 0 * stride;
  pixel[13] = -3 + 1 * stride;
  pixel[14] = -2 + 2 * stride;
  pixel[15] = -1 + 3 * stride;
  for(y = boundary ; y < ysize - boundary; y++)
  {
    cache_0 = im + boundary + y*stride;
    line_min = cache_0 - boundary;
    line_max = im + xsize - boundary + y * stride;

    cache_1 = cache_0 + pixel[5];
    cache_2 = cache_0 + pixel[14];

    for(; cache_0 < line_max;cache_0++, cache_1++, cache_2++)
    {
      cb = *cache_0 + barrier;
      c_b = *cache_0 - barrier;
      if(*cache_1 > cb)
        if(*cache_2 > cb)
          if(*(cache_0+3) > cb)
            if(*(cache_0 + pixel[0]) > cb)
              if(*(cache_0 + pixel[3]) > cb)
                if(*(cache_0 + pixel[6]) > cb)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[15]) > cb)
                      if(*(cache_0 + pixel[1]) > cb)
                        goto success;
                      else if(*(cache_0 + pixel[1]) < c_b)
                        continue;
                      else
                        if(*(cache_0 + pixel[10]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            if(*(cache_0 + pixel[7]) > cb)
                              if(*(cache_0 + pixel[9]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else if(*(cache_0 + pixel[15]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[8]) > cb)
                        if(*(cache_0 + pixel[7]) > cb)
                          if(*(cache_0 + pixel[1]) > cb)
                            goto success;
                          else if(*(cache_0 + pixel[1]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[10]) > cb)
                              goto success;
                            else
                              continue;
                        else
                          continue;
                      else
                        continue;
                  else if(*(cache_2+4) < c_b)
                    continue;
                  else
                    if(*(cache_1+-6) > cb)
                      if(*(cache_0 + pixel[9]) > cb)
                        if(*(cache_0 + pixel[10]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            if(*(cache_0 + pixel[7]) > cb)
                              goto success;
                            else if(*(cache_0 + pixel[7]) < c_b)
                              continue;
                            else
                              if(*(cache_0+-3) > cb)
                                goto success;
                              else
                                continue;
                          else if(*(cache_0 + pixel[8]) < c_b)
                            continue;
                          else
                            if(*(cache_0+-3) > cb)
                              if(*(cache_0 + pixel[1]) > cb)
                                if(*(cache_0 + pixel[13]) > cb)
                                  goto success;
                                else
                                  continue;
                              else
                                continue;
                            else
                              continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else if(*(cache_0 + pixel[6]) < c_b)
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*(cache_1+-6) > cb)
                        continue;
                      else if(*(cache_1+-6) < c_b)
                        if(*(cache_0 + pixel[15]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[15]) > cb)
                      if(*(cache_2+4) > cb)
                        if(*(cache_0 + pixel[1]) > cb)
                          goto success;
                        else if(*(cache_0 + pixel[1]) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[8]) > cb)
                            if(*(cache_1+-6) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                      else if(*(cache_2+4) < c_b)
                        continue;
                      else
                        if(*(cache_0 + pixel[9]) > cb)
                          if(*(cache_0+-3) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              if(*(cache_0 + pixel[10]) > cb)
                                if(*(cache_1+-6) > cb)
                                  goto success;
                                else
                                  continue;
                              else
                                continue;
                            else if(*(cache_0 + pixel[1]) < c_b)
                              continue;
                            else
                              if(*(cache_0 + pixel[8]) > cb)
                                if(*(cache_0 + pixel[10]) > cb)
                                  goto success;
                                else
                                  continue;
                              else
                                continue;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else
                    continue;
              else if(*(cache_0 + pixel[3]) < c_b)
                continue;
              else
                if(*(cache_0+-3) > cb)
                  if(*(cache_0 + pixel[10]) > cb)
                    if(*(cache_1+-6) > cb)
                      if(*(cache_0 + pixel[8]) > cb)
                        if(*(cache_0 + pixel[9]) > cb)
                          goto success;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          continue;
                        else
                          if(*(cache_2+4) > cb)
                            goto success;
                          else
                            continue;
                      else if(*(cache_0 + pixel[8]) < c_b)
                        if(*(cache_0 + pixel[7]) > cb || *(cache_0 + pixel[7]) < c_b)
                          continue;
                        else
                          goto success;
                      else
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[13]) > cb)
                            if(*(cache_0 + pixel[15]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else if(*(cache_2+4) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[9]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              if(*(cache_0 + pixel[13]) > cb)
                                if(*(cache_0 + pixel[15]) > cb)
                                  goto success;
                                else
                                  continue;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else if(*(cache_0 + pixel[0]) < c_b)
              if(*(cache_0 + pixel[7]) > cb)
                if(*(cache_0 + pixel[10]) > cb)
                  goto success;
                else
                  continue;
              else
                continue;
            else
              if(*(cache_0 + pixel[7]) > cb)
                if(*(cache_0 + pixel[10]) > cb)
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_0 + pixel[6]) > cb)
                      if(*(cache_0 + pixel[8]) > cb)
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[9]) > cb)
                            goto success;
                          else
                            continue;
                        else if(*(cache_2+4) < c_b)
                          continue;
                        else
                          if(*(cache_1+-6) > cb)
                            if(*(cache_0 + pixel[9]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                      else
                        continue;
                    else if(*(cache_0 + pixel[6]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[15]) > cb)
                        if(*(cache_0+-3) > cb)
                          if(*(cache_0 + pixel[9]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else if(*(cache_0 + pixel[3]) < c_b)
                    continue;
                  else
                    if(*(cache_0+-3) > cb)
                      if(*(cache_0 + pixel[8]) > cb)
                        if(*(cache_1+-6) > cb)
                          if(*(cache_0 + pixel[6]) > cb)
                            if(*(cache_0 + pixel[9]) > cb)
                              goto success;
                            else
                              continue;
                          else if(*(cache_0 + pixel[6]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[15]) > cb)
                              if(*(cache_0 + pixel[13]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else if(*(cache_0 + pixel[10]) < c_b)
                  continue;
                else
                  if(*(cache_0 + pixel[1]) > cb)
                    if(*(cache_0 + pixel[9]) > cb)
                      if(*(cache_0 + pixel[6]) > cb)
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[3]) > cb)
                            if(*(cache_0 + pixel[8]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else if(*(cache_0+3) < c_b)
            if(*(cache_0+-3) > cb)
              if(*(cache_0 + pixel[9]) > cb)
                if(*(cache_1+-6) > cb)
                  if(*(cache_0 + pixel[10]) > cb)
                    if(*(cache_0 + pixel[6]) > cb)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else
              continue;
          else
            if(*(cache_0+-3) > cb)
              if(*(cache_1+-6) > cb)
                if(*(cache_0 + pixel[7]) > cb)
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[9]) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            goto success;
                          else if(*(cache_0 + pixel[8]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[0]) > cb)
                              if(*(cache_0 + pixel[1]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                        else if(*(cache_0 + pixel[6]) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[15]) > cb)
                            if(*(cache_0 + pixel[8]) > cb)
                              goto success;
                            else if(*(cache_0 + pixel[8]) < c_b)
                              continue;
                            else
                              if(*(cache_0 + pixel[0]) > cb)
                                goto success;
                              else
                                continue;
                          else
                            continue;
                      else if(*(cache_0 + pixel[9]) < c_b)
                        continue;
                      else
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[0]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_0 + pixel[1]) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
                else if(*(cache_0 + pixel[7]) < c_b)
                  if(*(cache_0 + pixel[10]) > cb)
                    if(*(cache_2+4) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_0 + pixel[0]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[0]) > cb)
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else if(*(cache_2+4) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[9]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              if(*(cache_0 + pixel[15]) > cb)
                                goto success;
                              else
                                continue;
                            else if(*(cache_0 + pixel[1]) < c_b)
                              continue;
                            else
                              if(*(cache_0 + pixel[8]) > cb)
                                goto success;
                              else
                                continue;
                          else
                            continue;
                      else
                        continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            if(*(cache_0 + pixel[13]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else
                continue;
            else
              continue;
        else if(*cache_2 < c_b)
          if(*(cache_0+3) > cb)
            if(*(cache_0 + pixel[7]) > cb)
              if(*(cache_0 + pixel[1]) > cb)
                if(*(cache_0 + pixel[9]) > cb)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[6]) > cb)
                      if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_0 + pixel[8]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_2+4) < c_b)
                    continue;
                  else
                    if(*(cache_1+-6) > cb)
                      if(*(cache_0 + pixel[3]) > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                else if(*(cache_0 + pixel[9]) < c_b)
                  if(*(cache_0 + pixel[15]) > cb)
                    if(*(cache_2+4) > cb)
                      if(*(cache_0 + pixel[3]) > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[0]) > cb)
                    if(*(cache_0 + pixel[8]) > cb)
                      if(*(cache_2+4) > cb)
                        if(*(cache_0 + pixel[3]) > cb)
                          if(*(cache_0 + pixel[6]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else if(*(cache_0 + pixel[8]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[15]) > cb)
                        if(*(cache_2+4) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else if(*(cache_0 + pixel[1]) < c_b)
                if(*(cache_1+-6) > cb)
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[6]) > cb)
                        if(*(cache_0 + pixel[8]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_0 + pixel[3]) < c_b)
                    continue;
                  else
                    if(*(cache_0+-3) > cb)
                      if(*(cache_0 + pixel[10]) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else if(*(cache_1+-6) < c_b)
                  if(*(cache_0 + pixel[9]) > cb)
                    if(*(cache_0 + pixel[3]) > cb)
                      if(*(cache_2+4) > cb)
                        if(*(cache_0 + pixel[10]) > cb)
                          goto success;
                        else
                          continue;
                      else if(*(cache_2+4) < c_b)
                        if(*(cache_0 + pixel[10]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else if(*(cache_0 + pixel[3]) < c_b)
                      if(*(cache_0 + pixel[15]) < c_b)
                        if(*(cache_0+-3) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_2+4) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                  else if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0 + pixel[0]) < c_b)
                      goto success;
                    else
                      continue;
                  else
                    if(*(cache_2+4) < c_b)
                      if(*(cache_0 + pixel[10]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[15]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        if(*(cache_0 + pixel[3]) < c_b)
                          if(*(cache_0 + pixel[15]) < c_b)
                            if(*(cache_0 + pixel[0]) < c_b)
                              if(*(cache_0+-3) < c_b)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                else
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[8]) > cb)
                      goto success;
                    else
                      continue;
                  else
                    continue;
              else
                if(*(cache_0 + pixel[10]) > cb)
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_2+4) > cb)
                      if(*(cache_0 + pixel[6]) > cb)
                        if(*(cache_0 + pixel[9]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else if(*(cache_2+4) < c_b)
                      continue;
                    else
                      if(*(cache_1+-6) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                          if(*(cache_0 + pixel[9]) > cb)
                            if(*(cache_0 + pixel[8]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else if(*(cache_0 + pixel[3]) < c_b)
                    continue;
                  else
                    if(*(cache_0+-3) > cb)
                      if(*(cache_1+-6) > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                else
                  continue;
            else if(*(cache_0 + pixel[7]) < c_b)
              if(*(cache_1+-6) < c_b)
                if(*(cache_0 + pixel[15]) > cb)
                  continue;
                else if(*(cache_0 + pixel[15]) < c_b)
                  if(*(cache_0+-3) < c_b)
                    if(*(cache_0 + pixel[10]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0 + pixel[13]) < c_b)
                        if(*(cache_0 + pixel[9]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          if(*(cache_0 + pixel[8]) > cb)
                            continue;
                          else if(*(cache_0 + pixel[8]) < c_b)
                            goto success;
                          else
                            if(*(cache_0 + pixel[1]) < c_b)
                              goto success;
                            else
                              continue;
                        else
                          if(*(cache_2+4) < c_b)
                            goto success;
                          else
                            continue;
                      else
                        continue;
                    else
                      if(*(cache_0 + pixel[3]) < c_b)
                        goto success;
                      else
                        continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[6]) < c_b)
                    if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0+-3) < c_b)
                        if(*(cache_0 + pixel[8]) < c_b)
                          if(*(cache_0 + pixel[13]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
            else
              if(*(cache_0 + pixel[0]) < c_b)
                if(*(cache_0 + pixel[10]) > cb)
                  continue;
                else if(*(cache_0 + pixel[10]) < c_b)
                  if(*(cache_0 + pixel[9]) > cb)
                    continue;
                  else if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0+-3) < c_b)
                      if(*(cache_0 + pixel[1]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[1]) < c_b)
                        if(*(cache_1+-6) < c_b)
                          if(*(cache_0 + pixel[13]) < c_b)
                            if(*(cache_0 + pixel[15]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        if(*(cache_0 + pixel[8]) < c_b)
                          if(*(cache_0 + pixel[13]) < c_b)
                            if(*(cache_1+-6) < c_b)
                              if(*(cache_0 + pixel[15]) < c_b)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else
                    if(*(cache_2+4) < c_b)
                      if(*(cache_0+-3) < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          if(*(cache_1+-6) < c_b)
                            if(*(cache_0 + pixel[13]) < c_b)
                              if(*(cache_0 + pixel[15]) < c_b)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else
                  if(*(cache_0 + pixel[3]) < c_b)
                    if(*(cache_1+-6) < c_b)
                      if(*(cache_0+-3) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else if(*(cache_0+3) < c_b)
            if(*(cache_0+-3) > cb)
              if(*(cache_0 + pixel[13]) > cb)
                goto success;
              else
                continue;
            else if(*(cache_0+-3) < c_b)
              if(*(cache_0 + pixel[9]) > cb)
                if(*(cache_0 + pixel[13]) < c_b)
                  goto success;
                else
                  continue;
              else if(*(cache_0 + pixel[9]) < c_b)
                goto success;
              else
                if(*(cache_0 + pixel[6]) > cb || *(cache_0 + pixel[6]) < c_b)
                  continue;
                else
                  if(*(cache_2+4) < c_b)
                    goto success;
                  else
                    continue;
            else
              continue;
          else
            if(*(cache_1+-6) > cb)
              if(*(cache_0 + pixel[13]) > cb)
                if(*(cache_0 + pixel[9]) > cb)
                  if(*(cache_0 + pixel[7]) > cb)
                    if(*(cache_0+-3) > cb)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else if(*(cache_1+-6) < c_b)
              if(*(cache_0 + pixel[3]) > cb)
                if(*(cache_0 + pixel[8]) < c_b)
                  if(*(cache_0 + pixel[15]) > cb)
                    continue;
                  else if(*(cache_0 + pixel[15]) < c_b)
                    if(*(cache_0 + pixel[13]) < c_b)
                      if(*(cache_0 + pixel[0]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[0]) < c_b)
                        goto success;
                      else
                        if(*(cache_0 + pixel[7]) < c_b)
                          goto success;
                        else
                          continue;
                    else
                      continue;
                  else
                    if(*(cache_0 + pixel[6]) < c_b)
                      goto success;
                    else
                      continue;
                else
                  continue;
              else if(*(cache_0 + pixel[3]) < c_b)
                if(*(cache_2+4) > cb)
                  continue;
                else if(*(cache_2+4) < c_b)
                  if(*(cache_0 + pixel[0]) < c_b)
                    if(*(cache_0 + pixel[1]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[1]) < c_b)
                      if(*(cache_0 + pixel[15]) < c_b)
                        if(*(cache_0+-3) < c_b)
                          if(*(cache_0 + pixel[13]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      if(*(cache_0 + pixel[8]) < c_b)
                        goto success;
                      else
                        continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0 + pixel[1]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[1]) < c_b)
                      if(*(cache_0+-3) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      if(*(cache_0 + pixel[8]) < c_b)
                        if(*(cache_0 + pixel[0]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else
                if(*(cache_0 + pixel[1]) > cb)
                  continue;
                else if(*(cache_0 + pixel[1]) < c_b)
                  if(*(cache_0 + pixel[10]) < c_b)
                    if(*(cache_0+-3) < c_b)
                      if(*(cache_0 + pixel[9]) > cb)
                        if(*(cache_2+4) < c_b)
                          goto success;
                        else
                          continue;
                      else if(*(cache_0 + pixel[9]) < c_b)
                        if(*(cache_0 + pixel[15]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[15]) < c_b)
                          if(*(cache_0 + pixel[13]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          if(*(cache_0 + pixel[6]) < c_b)
                            goto success;
                          else
                            continue;
                      else
                        if(*(cache_2+4) < c_b)
                          if(*(cache_0 + pixel[15]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[7]) > cb)
                    continue;
                  else if(*(cache_0 + pixel[7]) < c_b)
                    if(*(cache_0 + pixel[15]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[15]) < c_b)
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0+-3) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            if(*(cache_0 + pixel[13]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      if(*(cache_0 + pixel[6]) < c_b)
                        goto success;
                      else
                        continue;
                  else
                    if(*(cache_0 + pixel[0]) < c_b)
                      if(*(cache_0 + pixel[8]) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
            else
              continue;
        else
          if(*(cache_0 + pixel[7]) > cb)
            if(*(cache_0 + pixel[3]) > cb)
              if(*(cache_0 + pixel[10]) > cb)
                if(*(cache_0+3) > cb)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[6]) > cb)
                      if(*(cache_0 + pixel[8]) > cb)
                        if(*(cache_0 + pixel[9]) > cb)
                          goto success;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[0]) > cb)
                            goto success;
                          else
                            continue;
                      else if(*(cache_0 + pixel[8]) < c_b)
                        continue;
                      else
                        if(*(cache_0 + pixel[15]) > cb)
                          if(*(cache_0 + pixel[0]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else if(*(cache_2+4) < c_b)
                    if(*(cache_1+-6) > cb)
                      goto success;
                    else
                      continue;
                  else
                    if(*(cache_1+-6) > cb)
                      if(*(cache_0 + pixel[6]) > cb)
                        if(*(cache_0 + pixel[8]) > cb)
                          if(*(cache_0 + pixel[9]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else if(*(cache_0+3) < c_b)
                  continue;
                else
                  if(*(cache_0+-3) > cb)
                    if(*(cache_0 + pixel[13]) > cb)
                      if(*(cache_1+-6) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            if(*(cache_0 + pixel[9]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else if(*(cache_0 + pixel[10]) < c_b)
                if(*(cache_0 + pixel[15]) > cb)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[6]) > cb)
                      if(*(cache_0+3) > cb)
                        if(*(cache_0 + pixel[0]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else if(*(cache_0 + pixel[15]) < c_b)
                  continue;
                else
                  if(*(cache_0 + pixel[8]) > cb)
                    if(*(cache_0 + pixel[0]) > cb)
                      goto success;
                    else if(*(cache_0 + pixel[0]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[9]) > cb)
                        if(*(cache_0 + pixel[1]) > cb)
                          if(*(cache_0 + pixel[6]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else
                if(*(cache_0 + pixel[1]) > cb)
                  if(*(cache_0 + pixel[9]) > cb)
                    if(*(cache_0 + pixel[6]) > cb)
                      if(*(cache_0+3) > cb)
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            goto success;
                          else if(*(cache_0 + pixel[8]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[15]) > cb)
                              goto success;
                            else
                              continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0 + pixel[0]) > cb)
                      goto success;
                    else
                      continue;
                  else
                    if(*(cache_0 + pixel[0]) > cb)
                      if(*(cache_0+3) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            if(*(cache_2+4) > cb)
                              goto success;
                            else
                              continue;
                          else if(*(cache_0 + pixel[15]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[8]) > cb)
                              if(*(cache_2+4) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else
                  continue;
            else if(*(cache_0 + pixel[3]) < c_b)
              if(*(cache_0 + pixel[13]) > cb)
                if(*(cache_1+-6) > cb)
                  if(*(cache_0 + pixel[9]) > cb)
                    if(*(cache_0+-3) > cb)
                      if(*(cache_0 + pixel[6]) > cb)
                        if(*(cache_0 + pixel[8]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else if(*(cache_0 + pixel[13]) < c_b)
                continue;
              else
                if(*(cache_0+3) > cb)
                  if(*(cache_0+-3) > cb)
                    if(*(cache_0 + pixel[10]) > cb)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else
              if(*(cache_0+-3) > cb)
                if(*(cache_0 + pixel[13]) > cb)
                  if(*(cache_1+-6) > cb)
                    if(*(cache_0 + pixel[9]) > cb)
                      if(*(cache_0 + pixel[6]) > cb)
                        if(*(cache_0 + pixel[10]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else if(*(cache_0 + pixel[13]) < c_b)
                  if(*(cache_0 + pixel[0]) > cb)
                    goto success;
                  else
                    continue;
                else
                  if(*(cache_0+3) > cb)
                    if(*(cache_0 + pixel[9]) > cb)
                      if(*(cache_1+-6) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                          if(*(cache_0 + pixel[10]) > cb)
                            if(*(cache_0 + pixel[8]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else
            continue;
      else if(*cache_1 < c_b)
        if(*(cache_0 + pixel[15]) > cb)
          if(*(cache_1+-6) > cb)
            if(*(cache_2+4) > cb)
              if(*(cache_0+-3) > cb)
                if(*(cache_0 + pixel[10]) > cb)
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*cache_2 > cb)
                        goto success;
                      else
                        continue;
                    else if(*(cache_0 + pixel[1]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[7]) > cb)
                        goto success;
                      else
                        continue;
                  else
                    continue;
                else if(*(cache_0 + pixel[10]) < c_b)
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_0 + pixel[13]) > cb)
                      if(*cache_2 > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*cache_2 > cb)
                        if(*(cache_0 + pixel[0]) > cb)
                          if(*(cache_0 + pixel[13]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
            else if(*(cache_2+4) < c_b)
              if(*(cache_0 + pixel[7]) > cb)
                if(*(cache_0+-3) > cb)
                  if(*cache_2 > cb)
                    if(*(cache_0 + pixel[13]) > cb)
                      if(*(cache_0 + pixel[9]) > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else if(*(cache_0 + pixel[7]) < c_b)
                if(*(cache_0 + pixel[9]) > cb)
                  if(*(cache_0 + pixel[1]) > cb)
                    if(*(cache_0+-3) > cb)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else if(*(cache_0 + pixel[9]) < c_b)
                  if(*(cache_0 + pixel[10]) > cb)
                    continue;
                  else if(*(cache_0 + pixel[10]) < c_b)
                    if(*(cache_0 + pixel[3]) < c_b)
                      if(*(cache_0+3) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    if(*(cache_0 + pixel[1]) < c_b)
                      if(*(cache_0 + pixel[3]) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                else
                  if(*(cache_0 + pixel[0]) < c_b)
                    goto success;
                  else
                    continue;
              else
                if(*(cache_0 + pixel[0]) > cb)
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[9]) > cb)
                      if(*cache_2 > cb)
                        if(*(cache_0 + pixel[1]) > cb)
                          if(*(cache_0 + pixel[10]) > cb)
                            goto success;
                          else
                            continue;
                        else if(*(cache_0 + pixel[1]) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[8]) > cb)
                            if(*(cache_0+-3) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else
              if(*(cache_0 + pixel[9]) > cb)
                if(*(cache_0+-3) > cb)
                  if(*(cache_0 + pixel[1]) > cb)
                    if(*cache_2 > cb)
                      if(*(cache_0 + pixel[10]) > cb)
                        if(*(cache_0 + pixel[13]) > cb)
                          if(*(cache_0 + pixel[0]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_0 + pixel[1]) < c_b)
                    continue;
                  else
                    if(*(cache_0 + pixel[7]) > cb)
                      if(*(cache_0 + pixel[10]) > cb)
                        if(*(cache_0 + pixel[13]) > cb)
                          if(*cache_2 > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else if(*(cache_0 + pixel[7]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[0]) > cb)
                        if(*(cache_0 + pixel[8]) > cb)
                          if(*(cache_0 + pixel[6]) < c_b)
                            if(*(cache_0 + pixel[10]) > cb)
                              if(*(cache_0 + pixel[13]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            goto success;
                        else
                          continue;
                      else
                        continue;
                else
                  continue;
              else
                continue;
          else if(*(cache_1+-6) < c_b)
            if(*(cache_0 + pixel[3]) > cb)
              if(*(cache_0 + pixel[13]) > cb)
                if(*(cache_0+-3) > cb)
                  if(*(cache_0+3) > cb)
                    goto success;
                  else
                    continue;
                else if(*(cache_0+-3) < c_b)
                  if(*(cache_0+3) < c_b)
                    if(*(cache_0 + pixel[6]) < c_b)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else if(*(cache_0 + pixel[13]) < c_b)
                if(*(cache_0 + pixel[7]) < c_b)
                  if(*(cache_0 + pixel[6]) < c_b)
                    if(*(cache_0 + pixel[8]) < c_b)
                      if(*(cache_0+-3) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                if(*(cache_0+3) < c_b)
                  if(*(cache_0+-3) < c_b)
                    if(*(cache_0 + pixel[7]) < c_b)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else if(*(cache_0 + pixel[3]) < c_b)
              if(*(cache_0 + pixel[8]) < c_b)
                if(*(cache_0 + pixel[9]) < c_b)
                  if(*(cache_0 + pixel[7]) < c_b)
                    if(*(cache_0+3) > cb)
                      continue;
                    else if(*(cache_0+3) < c_b)
                      if(*(cache_0 + pixel[10]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[6]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        if(*(cache_0 + pixel[1]) < c_b)
                          goto success;
                        else
                          continue;
                    else
                      if(*(cache_0 + pixel[13]) < c_b)
                        goto success;
                      else
                        continue;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else
              if(*(cache_0+-3) < c_b)
                if(*(cache_0+3) > cb)
                  continue;
                else if(*(cache_0+3) < c_b)
                  if(*(cache_0 + pixel[6]) < c_b)
                    if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0 + pixel[9]) < c_b)
                        if(*(cache_0 + pixel[7]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[13]) < c_b)
                    if(*(cache_0 + pixel[7]) < c_b)
                      if(*(cache_0 + pixel[6]) < c_b)
                        if(*(cache_0 + pixel[10]) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            if(*(cache_0 + pixel[9]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else
            if(*(cache_2+4) > cb)
              if(*(cache_0+3) > cb)
                if(*(cache_0+-3) > cb)
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*(cache_0 + pixel[3]) > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else if(*(cache_2+4) < c_b)
              if(*(cache_0 + pixel[10]) > cb)
                continue;
              else if(*(cache_0 + pixel[10]) < c_b)
                if(*(cache_0+3) < c_b)
                  if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0 + pixel[3]) < c_b)
                      if(*(cache_0 + pixel[7]) < c_b)
                        if(*(cache_0 + pixel[6]) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                if(*(cache_0 + pixel[1]) < c_b)
                  if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0 + pixel[3]) < c_b)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else
              continue;
        else if(*(cache_0 + pixel[15]) < c_b)
          if(*(cache_0+3) > cb)
            if(*(cache_0+-3) < c_b)
              if(*(cache_1+-6) < c_b)
                if(*(cache_0 + pixel[13]) < c_b)
                  if(*(cache_0 + pixel[7]) > cb)
                    continue;
                  else if(*(cache_0 + pixel[7]) < c_b)
                    goto success;
                  else
                    if(*(cache_0 + pixel[8]) < c_b)
                      if(*(cache_0 + pixel[0]) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                else
                  continue;
              else
                continue;
            else
              continue;
          else if(*(cache_0+3) < c_b)
            if(*(cache_0 + pixel[6]) > cb)
              if(*(cache_0 + pixel[13]) > cb)
                if(*cache_2 > cb)
                  if(*(cache_0 + pixel[10]) > cb)
                    goto success;
                  else
                    continue;
                else
                  continue;
              else if(*(cache_0 + pixel[13]) < c_b)
                if(*(cache_0 + pixel[0]) < c_b)
                  if(*(cache_2+4) < c_b)
                    if(*cache_2 < c_b)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else if(*(cache_0 + pixel[6]) < c_b)
              if(*(cache_0 + pixel[3]) > cb)
                if(*(cache_0+-3) < c_b)
                  if(*(cache_0 + pixel[1]) < c_b)
                    continue;
                  else
                    goto success;
                else
                  continue;
              else if(*(cache_0 + pixel[3]) < c_b)
                if(*(cache_0 + pixel[7]) > cb)
                  if(*cache_2 < c_b)
                    goto success;
                  else
                    continue;
                else if(*(cache_0 + pixel[7]) < c_b)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[10]) < c_b)
                      goto success;
                    else
                      continue;
                  else if(*(cache_2+4) < c_b)
                    if(*(cache_0 + pixel[1]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[1]) < c_b)
                      if(*(cache_0 + pixel[0]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[0]) < c_b)
                        goto success;
                      else
                        if(*(cache_0 + pixel[9]) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                    else
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[8]) < c_b)
                          if(*(cache_0 + pixel[9]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    if(*(cache_1+-6) < c_b)
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[8]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else
                  if(*cache_2 < c_b)
                    if(*(cache_2+4) < c_b)
                      if(*(cache_0 + pixel[0]) < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                if(*(cache_0+-3) < c_b)
                  if(*(cache_1+-6) < c_b)
                    if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0 + pixel[8]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[8]) < c_b)
                        if(*(cache_0 + pixel[9]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          if(*(cache_0 + pixel[7]) > cb)
                            continue;
                          else if(*(cache_0 + pixel[7]) < c_b)
                            goto success;
                          else
                            if(*(cache_0 + pixel[13]) < c_b)
                              goto success;
                            else
                              continue;
                        else
                          if(*(cache_2+4) < c_b)
                            goto success;
                          else
                            continue;
                      else
                        if(*(cache_0 + pixel[13]) < c_b)
                          if(*(cache_0 + pixel[0]) < c_b)
                            if(*(cache_0 + pixel[7]) > cb || *(cache_0 + pixel[7]) < c_b)
                              continue;
                            else
                              goto success;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else
              if(*(cache_0 + pixel[13]) < c_b)
                if(*(cache_2+4) > cb)
                  continue;
                else if(*(cache_2+4) < c_b)
                  if(*cache_2 < c_b)
                    if(*(cache_0 + pixel[3]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[3]) < c_b)
                      if(*(cache_0 + pixel[0]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[0]) < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        if(*(cache_0 + pixel[7]) < c_b)
                          if(*(cache_1+-6) < c_b)
                            if(*(cache_0 + pixel[8]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else
                      if(*(cache_0+-3) < c_b)
                        if(*(cache_0 + pixel[10]) < c_b)
                          if(*(cache_0 + pixel[1]) > cb)
                            continue;
                          else if(*(cache_0 + pixel[1]) < c_b)
                            goto success;
                          else
                            if(*(cache_0 + pixel[7]) < c_b)
                              goto success;
                            else
                              continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
                else
                  if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_1+-6) < c_b)
                      if(*(cache_0 + pixel[0]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[0]) < c_b)
                        if(*cache_2 < c_b)
                          if(*(cache_0 + pixel[10]) < c_b)
                            if(*(cache_0+-3) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        if(*(cache_0 + pixel[7]) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            if(*(cache_0 + pixel[1]) > cb || *(cache_0 + pixel[1]) < c_b)
                              continue;
                            else
                              goto success;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else
            if(*(cache_0+-3) < c_b)
              if(*(cache_0 + pixel[13]) < c_b)
                if(*(cache_1+-6) < c_b)
                  if(*(cache_0 + pixel[9]) > cb)
                    if(*(cache_0 + pixel[3]) < c_b)
                      if(*(cache_2+4) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_0 + pixel[9]) < c_b)
                    if(*(cache_0 + pixel[10]) > cb)
                      continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0 + pixel[7]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[7]) < c_b)
                        if(*cache_2 > cb || *cache_2 < c_b)
                          goto success;
                        else
                          if(*(cache_0 + pixel[6]) < c_b)
                            if(*(cache_0 + pixel[8]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                      else
                        if(*(cache_0 + pixel[1]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[1]) < c_b)
                          if(*cache_2 < c_b)
                            if(*(cache_0 + pixel[0]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          if(*(cache_0 + pixel[0]) < c_b)
                            if(*(cache_0 + pixel[8]) < c_b)
                              if(*cache_2 < c_b)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                    else
                      if(*(cache_0 + pixel[3]) < c_b)
                        goto success;
                      else
                        continue;
                  else
                    if(*(cache_2+4) < c_b)
                      if(*(cache_0 + pixel[1]) < c_b)
                        if(*(cache_0 + pixel[10]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[10]) < c_b)
                          if(*cache_2 < c_b)
                            if(*(cache_0 + pixel[0]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          if(*(cache_0 + pixel[3]) < c_b)
                            if(*(cache_0 + pixel[0]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                      else
                        continue;
                    else
                      continue;
                else
                  continue;
              else
                continue;
            else
              continue;
        else
          if(*(cache_0 + pixel[8]) > cb)
            if(*(cache_0 + pixel[6]) > cb)
              if(*cache_2 > cb)
                if(*(cache_1+-6) > cb)
                  if(*(cache_0 + pixel[10]) > cb)
                    goto success;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else
              continue;
          else if(*(cache_0 + pixel[8]) < c_b)
            if(*(cache_0 + pixel[3]) > cb)
              if(*(cache_0 + pixel[13]) > cb)
                continue;
              else if(*(cache_0 + pixel[13]) < c_b)
                if(*(cache_0+-3) < c_b)
                  if(*(cache_0 + pixel[7]) < c_b)
                    if(*(cache_1+-6) < c_b)
                      if(*(cache_0 + pixel[6]) < c_b)
                        if(*(cache_0 + pixel[10]) < c_b)
                          if(*(cache_0 + pixel[9]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                if(*(cache_0+3) < c_b)
                  if(*(cache_0+-3) < c_b)
                    if(*(cache_0 + pixel[10]) < c_b)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else if(*(cache_0 + pixel[3]) < c_b)
              if(*(cache_2+4) > cb)
                if(*(cache_1+-6) < c_b)
                  if(*(cache_0 + pixel[7]) < c_b)
                    goto success;
                  else
                    continue;
                else
                  continue;
              else if(*(cache_2+4) < c_b)
                if(*(cache_0 + pixel[6]) < c_b)
                  if(*(cache_0+3) > cb)
                    continue;
                  else if(*(cache_0+3) < c_b)
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[0]) > cb)
                        continue;
                      else if(*(cache_0 + pixel[0]) < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          goto success;
                        else
                          continue;
                      else
                        if(*(cache_0 + pixel[9]) < c_b)
                          if(*(cache_0 + pixel[1]) < c_b)
                            if(*(cache_0 + pixel[7]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0 + pixel[7]) < c_b)
                        if(*(cache_0 + pixel[9]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          goto success;
                        else
                          if(*(cache_0 + pixel[0]) < c_b)
                            goto success;
                          else
                            continue;
                      else
                        continue;
                    else
                      if(*(cache_0 + pixel[1]) < c_b)
                        if(*(cache_0 + pixel[9]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          if(*(cache_0 + pixel[7]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          if(*(cache_0 + pixel[0]) < c_b)
                            if(*(cache_0 + pixel[7]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                      else
                        continue;
                  else
                    if(*(cache_0+-3) < c_b)
                      if(*(cache_0 + pixel[13]) < c_b)
                        if(*(cache_1+-6) < c_b)
                          if(*(cache_0 + pixel[7]) < c_b)
                            if(*(cache_0 + pixel[10]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else
                  continue;
              else
                if(*(cache_1+-6) < c_b)
                  if(*(cache_0+3) > cb)
                    continue;
                  else if(*(cache_0+3) < c_b)
                    if(*(cache_0 + pixel[6]) < c_b)
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[7]) < c_b)
                          if(*(cache_0 + pixel[9]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    if(*(cache_0+-3) < c_b)
                      if(*(cache_0 + pixel[13]) < c_b)
                        if(*(cache_0 + pixel[6]) < c_b)
                          if(*(cache_0 + pixel[7]) < c_b)
                            if(*(cache_0 + pixel[10]) < c_b)
                              if(*(cache_0 + pixel[9]) < c_b)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else
                  continue;
            else
              if(*(cache_0+-3) < c_b)
                if(*(cache_0 + pixel[13]) > cb)
                  if(*(cache_0+3) < c_b)
                    goto success;
                  else
                    continue;
                else if(*(cache_0 + pixel[13]) < c_b)
                  if(*(cache_1+-6) < c_b)
                    if(*(cache_0 + pixel[7]) < c_b)
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[6]) < c_b)
                          if(*(cache_0 + pixel[9]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0+3) < c_b)
                    if(*(cache_0 + pixel[10]) < c_b)
                      if(*(cache_0 + pixel[6]) < c_b)
                        if(*(cache_1+-6) < c_b)
                          if(*(cache_0 + pixel[7]) < c_b)
                            if(*(cache_0 + pixel[9]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else
            continue;
      else
        if(*(cache_0+-3) > cb)
          if(*cache_2 > cb)
            if(*(cache_0 + pixel[7]) > cb)
              if(*(cache_1+-6) > cb)
                if(*(cache_0 + pixel[6]) > cb)
                  if(*(cache_0 + pixel[13]) > cb)
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[9]) > cb)
                        if(*(cache_0 + pixel[8]) > cb)
                          goto success;
                        else if(*(cache_0 + pixel[8]) < c_b)
                          continue;
                        else
                          if(*(cache_0 + pixel[0]) > cb)
                            goto success;
                          else
                            continue;
                      else if(*(cache_0 + pixel[9]) < c_b)
                        continue;
                      else
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[0]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_0 + pixel[0]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
                else if(*(cache_0 + pixel[6]) < c_b)
                  continue;
                else
                  if(*(cache_0 + pixel[15]) > cb)
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_0 + pixel[9]) > cb)
                          if(*(cache_0 + pixel[8]) > cb)
                            goto success;
                          else if(*(cache_0 + pixel[8]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[1]) > cb)
                              goto success;
                            else
                              continue;
                        else if(*(cache_0 + pixel[9]) < c_b)
                          continue;
                        else
                          if(*(cache_2+4) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              if(*(cache_0 + pixel[0]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                      else
                        continue;
                    else if(*(cache_0 + pixel[10]) < c_b)
                      continue;
                    else
                      if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_0 + pixel[1]) > cb)
                          if(*(cache_2+4) > cb)
                            if(*(cache_0 + pixel[13]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else if(*(cache_1+-6) < c_b)
                continue;
              else
                if(*(cache_0+3) > cb)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*(cache_0 + pixel[0]) > cb)
                        if(*(cache_0 + pixel[3]) > cb)
                          if(*(cache_0 + pixel[13]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else if(*(cache_0 + pixel[7]) < c_b)
              if(*(cache_2+4) > cb)
                if(*(cache_1+-6) > cb)
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_0 + pixel[15]) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_0 + pixel[1]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_0 + pixel[3]) < c_b)
                    continue;
                  else
                    if(*(cache_0 + pixel[10]) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_0 + pixel[0]) > cb)
                          if(*(cache_0 + pixel[1]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                else if(*(cache_1+-6) < c_b)
                  if(*(cache_0+3) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*(cache_0 + pixel[0]) > cb)
                        if(*(cache_0 + pixel[3]) > cb)
                          goto success;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0+3) > cb)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_0 + pixel[3]) > cb)
                          if(*(cache_0 + pixel[0]) > cb)
                            if(*(cache_0 + pixel[15]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else if(*(cache_2+4) < c_b)
                continue;
              else
                if(*(cache_0 + pixel[9]) > cb)
                  if(*(cache_0 + pixel[0]) > cb)
                    if(*(cache_1+-6) > cb)
                      goto success;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
            else
              if(*(cache_0 + pixel[0]) > cb)
                if(*(cache_0 + pixel[10]) > cb)
                  if(*(cache_2+4) > cb)
                    if(*(cache_0 + pixel[13]) > cb)
                      if(*(cache_1+-6) > cb)
                        if(*(cache_0 + pixel[15]) > cb)
                          if(*(cache_0 + pixel[1]) > cb)
                            goto success;
                          else if(*(cache_0 + pixel[1]) < c_b)
                            continue;
                          else
                            if(*(cache_0 + pixel[8]) > cb)
                              goto success;
                            else
                              continue;
                        else
                          continue;
                      else if(*(cache_1+-6) < c_b)
                        continue;
                      else
                        if(*(cache_0+3) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                  else if(*(cache_2+4) < c_b)
                    if(*(cache_0 + pixel[1]) > cb)
                      if(*(cache_0 + pixel[3]) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                  else
                    if(*(cache_0 + pixel[9]) > cb)
                      if(*(cache_0 + pixel[1]) > cb)
                        if(*(cache_0 + pixel[13]) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            if(*(cache_1+-6) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else if(*(cache_0 + pixel[1]) < c_b)
                        continue;
                      else
                        if(*(cache_0 + pixel[8]) > cb)
                          if(*(cache_1+-6) > cb)
                            if(*(cache_0 + pixel[13]) > cb)
                              if(*(cache_0 + pixel[15]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                    else
                      continue;
                else if(*(cache_0 + pixel[10]) < c_b)
                  if(*(cache_0+3) > cb)
                    if(*(cache_0 + pixel[13]) > cb)
                      if(*(cache_2+4) > cb)
                        if(*(cache_0 + pixel[3]) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else if(*(cache_0+3) < c_b)
                    continue;
                  else
                    if(*(cache_1+-6) > cb)
                      if(*(cache_0 + pixel[3]) > cb)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                else
                  if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_1+-6) > cb)
                      if(*(cache_0 + pixel[13]) > cb)
                        if(*(cache_2+4) > cb)
                          if(*(cache_0 + pixel[15]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else if(*(cache_1+-6) < c_b)
                      if(*(cache_0+3) > cb)
                        goto success;
                      else
                        continue;
                    else
                      if(*(cache_0+3) > cb)
                        if(*(cache_0 + pixel[13]) > cb)
                          if(*(cache_0 + pixel[1]) > cb)
                            if(*(cache_2+4) > cb)
                              if(*(cache_0 + pixel[15]) > cb)
                                goto success;
                              else
                                continue;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else
                continue;
          else
            continue;
        else if(*(cache_0+-3) < c_b)
          if(*(cache_0 + pixel[15]) > cb)
            if(*cache_2 < c_b)
              if(*(cache_0 + pixel[6]) < c_b)
                if(*(cache_0 + pixel[10]) < c_b)
                  if(*(cache_0 + pixel[7]) < c_b)
                    goto success;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else
              continue;
          else if(*(cache_0 + pixel[15]) < c_b)
            if(*(cache_0 + pixel[10]) > cb)
              if(*(cache_0+3) > cb)
                continue;
              else if(*(cache_0+3) < c_b)
                if(*(cache_0 + pixel[3]) < c_b)
                  if(*(cache_0 + pixel[13]) < c_b)
                    goto success;
                  else
                    continue;
                else
                  continue;
              else
                if(*(cache_1+-6) < c_b)
                  if(*(cache_0 + pixel[3]) < c_b)
                    goto success;
                  else
                    continue;
                else
                  continue;
            else if(*(cache_0 + pixel[10]) < c_b)
              if(*cache_2 < c_b)
                if(*(cache_0 + pixel[9]) > cb)
                  if(*(cache_2+4) < c_b)
                    goto success;
                  else
                    continue;
                else if(*(cache_0 + pixel[9]) < c_b)
                  if(*(cache_1+-6) > cb)
                    continue;
                  else if(*(cache_1+-6) < c_b)
                    if(*(cache_0 + pixel[13]) < c_b)
                      if(*(cache_0 + pixel[1]) > cb)
                        if(*(cache_0 + pixel[7]) < c_b)
                          goto success;
                        else
                          continue;
                      else if(*(cache_0 + pixel[1]) < c_b)
                        if(*(cache_0 + pixel[0]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[0]) < c_b)
                          goto success;
                        else
                          if(*(cache_0 + pixel[7]) < c_b)
                            goto success;
                          else
                            continue;
                      else
                        if(*(cache_0 + pixel[7]) > cb)
                          continue;
                        else if(*(cache_0 + pixel[7]) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          if(*(cache_0 + pixel[0]) < c_b)
                            if(*(cache_0 + pixel[8]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                    else
                      continue;
                  else
                    if(*(cache_0+3) < c_b)
                      if(*(cache_0 + pixel[3]) < c_b)
                        goto success;
                      else
                        continue;
                    else
                      continue;
                else
                  if(*(cache_2+4) < c_b)
                    if(*(cache_1+-6) > cb)
                      continue;
                    else if(*(cache_1+-6) < c_b)
                      if(*(cache_0 + pixel[13]) < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          if(*(cache_0 + pixel[0]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      if(*(cache_0+3) < c_b)
                        if(*(cache_0 + pixel[3]) < c_b)
                          if(*(cache_0 + pixel[0]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                  else
                    continue;
              else
                continue;
            else
              if(*(cache_0 + pixel[3]) < c_b)
                if(*(cache_1+-6) > cb)
                  continue;
                else if(*(cache_1+-6) < c_b)
                  if(*(cache_2+4) < c_b)
                    if(*(cache_0 + pixel[13]) < c_b)
                      if(*cache_2 < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          if(*(cache_0 + pixel[0]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  if(*(cache_0+3) < c_b)
                    if(*(cache_2+4) < c_b)
                      if(*cache_2 < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                          if(*(cache_0 + pixel[13]) < c_b)
                            if(*(cache_0 + pixel[0]) < c_b)
                              goto success;
                            else
                              continue;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
              else
                continue;
          else
            if(*(cache_0 + pixel[6]) < c_b)
              if(*cache_2 < c_b)
                if(*(cache_0 + pixel[7]) < c_b)
                  if(*(cache_1+-6) < c_b)
                    if(*(cache_0 + pixel[13]) < c_b)
                      if(*(cache_0 + pixel[10]) < c_b)
                        if(*(cache_0 + pixel[9]) < c_b)
                          if(*(cache_0 + pixel[8]) < c_b)
                            goto success;
                          else
                            continue;
                        else
                          continue;
                      else
                        continue;
                    else
                      continue;
                  else
                    continue;
                else
                  continue;
              else
                continue;
            else
              continue;
        else
          continue;
success:
      if(total >= rsize) {
        rsize *=2;
        ret=(xy*)realloc(ret, rsize*sizeof(xy));
      }
      ret[total].x = cache_0-line_min;
      ret[total++].y = y;
    }
  }
  *num = total;
  return ret;
}

int corner_score(const byte*  imp, const int *pointer_dir, int barrier) {
  /*The score for a positive feature is sum of the difference between the pixels
    and the barrier if the difference is positive. Negative is similar.
    The score is the max of those two.

    B = {x | x = points on the Bresenham circle around c}
    Sp = { I(x) - t | x E B , I(x) - t > 0 }
    Sn = { t - I(x) | x E B, t - I(x) > 0}

    Score = max sum(Sp), sum(Sn)*/

  int cb = *imp + barrier;
  int c_b = *imp - barrier;
  int sp=0, sn = 0;

  int i=0;

  for(i=0; i<16; i++)
  {
    int p = imp[pointer_dir[i]];

    if(p > cb)
      sp += p-cb;
    else if(p < c_b)
      sn += c_b-p;
  }

  if(sp > sn)
    return sp;
  else
    return sn;
}

xy* fast_nonmax(const byte* im, int xsize, int ysize, int stride,
                xy* corners, int numcorners, int barrier, int* numnx) {
  /*Create a list of integer pointer offstes, corresponding to the */
  /*direction offsets in dir[]*/
  int  pointer_dir[16];
  int* row_start = (int*) malloc(ysize * sizeof(int));
  int* scores    = (int*) malloc(numcorners * sizeof(int));
  xy*  nonmax_corners=(xy*)malloc(numcorners* sizeof(xy));
  int num_nonmax=0;
  int prev_row = -1;
  int i, j;
  int point_above = 0;
  int point_below = 0;
  (void) xsize;

  pointer_dir[0] = 0 + 3 * stride;
  pointer_dir[1] = 1 + 3 * stride;
  pointer_dir[2] = 2 + 2 * stride;
  pointer_dir[3] = 3 + 1 * stride;
  pointer_dir[4] = 3 + 0 * stride;
  pointer_dir[5] = 3 + -1 * stride;
  pointer_dir[6] = 2 + -2 * stride;
  pointer_dir[7] = 1 + -3 * stride;
  pointer_dir[8] = 0 + -3 * stride;
  pointer_dir[9] = -1 + -3 * stride;
  pointer_dir[10] = -2 + -2 * stride;
  pointer_dir[11] = -3 + -1 * stride;
  pointer_dir[12] = -3 + 0 * stride;
  pointer_dir[13] = -3 + 1 * stride;
  pointer_dir[14] = -2 + 2 * stride;
  pointer_dir[15] = -1 + 3 * stride;

  if (numcorners < 5) {
    free(row_start);
    free(scores);
    free(nonmax_corners);
    return 0;
  }

  /*xsize ysize numcorners corners*/
  /*Compute the score for each detected corner, and find where each row begins*/
  /* (the corners are output in raster scan order). A beginning of -1 signifies*/
  /* that there are no corners on that row.*/


  for (i = 0; i < ysize; i++)
    row_start[i] = -1;


  for (i = 0; i< numcorners; i++) {
    if (corners[i].y != prev_row) {
      row_start[corners[i].y] = i;
      prev_row = corners[i].y;
    }
    scores[i] = corner_score(im + corners[i].x + corners[i].y * stride,
                             pointer_dir, barrier);
  }
  /*Point above points (roughly) to the pixel above the one of interest, if there*/
  /*is a feature there.*/
  for (i = 0; i < numcorners; i++) {
    int score = scores[i];
    xy pos = corners[i];

    //Check left
    if (i > 0)
      if (corners[i-1].x == pos.x-1 && corners[i-1].y == pos.y &&
          scores[i-1] > score)
        continue;

    //Check right
    if (i < (numcorners - 1))
      if (corners[i+1].x == pos.x+1 && corners[i+1].y == pos.y && scores[i+1] > score)
        continue;

    //Check above (if there is a valid row above)
    if (pos.y != 0 && row_start[pos.y - 1] != -1) {
      //Make sure that current point_above is one
      //row above.
      if (corners[point_above].y < pos.y - 1)
        point_above = row_start[pos.y-1];

      //Make point_above point to the first of the pixels above the current point,
      //if it exists.
      for (; corners[point_above].y < pos.y &&
           corners[point_above].x < pos.x - 1; point_above++) {}


      for (j = point_above;
           corners[j].y < pos.y && corners[j].x <= pos.x + 1; j++) {
        int x = corners[j].x;
        if ((x == pos.x - 1 || x ==pos.x || x == pos.x+1) && (scores[j] > score))
          goto cont;
      }
    }

    //Check below (if there is anything below)
    if (pos.y != ysize-1 && row_start[pos.y + 1] != -1 &&
        point_below < numcorners) {  // Nothing below
      if (corners[point_below].y < pos.y + 1)
        point_below = row_start[pos.y+1];

      // Make point below point to one of the pixels belowthe current point,
      // if it exists.
      for (; point_below < numcorners && corners[point_below].y == pos.y+1 &&
             corners[point_below].x < pos.x - 1; point_below++) {}

      for(j = point_below;
          j < numcorners && corners[j].y == pos.y+1 &&
          corners[j].x <= pos.x + 1; j++) {
        int x = corners[j].x;
        if ((x == pos.x - 1 || x ==pos.x || x == pos.x+1) && (scores[j] >score))
          goto cont;
      }
    }

    nonmax_corners[num_nonmax].x = corners[i].x;
    nonmax_corners[num_nonmax].y = corners[i].y;
    num_nonmax++;

cont:
    ;
  }
  *numnx = num_nonmax;
  free(row_start);
  free(scores);
  return nonmax_corners;
}

//////////////////////////////////////////////////////////////////////////////
// Nonmax suppression wrapper

#define NONMAX_BARRIER  5
int NonmaxSuppression(unsigned char *frmbuf, int width, int height, int stride,
                      int *frm_corners, int num_frm_corners) {
  int num_frm_corners_nonmax;
  xy *frm_corners_nonmax_xy;
  frm_corners_nonmax_xy = fast_nonmax(frmbuf, width, height, stride,
                                      (xy *)frm_corners, num_frm_corners,
                                      NONMAX_BARRIER, &num_frm_corners_nonmax);
  if (frm_corners_nonmax_xy &&
      num_frm_corners_nonmax <= num_frm_corners) {
    memcpy(frm_corners, frm_corners_nonmax_xy,
           sizeof(xy) * num_frm_corners_nonmax);
    free(frm_corners_nonmax_xy);
    return num_frm_corners_nonmax;
  } else {
    return num_frm_corners;
  }
}

//////////////////////////////////////////////////////////////////////////////
// Fast_9 wrapper

#define FAST_BARRIER    40
int FastCornerDetect(unsigned char *buf, int width, int height, int stride,
                     int *points, int max_points) {
  int num_points;
  xy *frm_corners_xy;
  frm_corners_xy = fast_corner_detect_9(buf, width, height, stride,
                                        FAST_BARRIER, &num_points);
  num_points =
      (num_points <= max_points ? num_points : max_points);
  if (num_points > 0 && frm_corners_xy) {
    memcpy(points, frm_corners_xy, sizeof(xy) * num_points);
    free(frm_corners_xy);
    return NonmaxSuppression(buf, width, height, stride, points, num_points);
  } else {
    return 0;
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Harris Corner detector

#define HARRIS_THRESH   50

#define GAUSS_FILT_SZ 7
#define GAUSS_FILT_SZ_BY2 ((GAUSS_FILT_SZ - 1) / 2)

double GaussianKernel[GAUSS_FILT_SZ_BY2 + 1][GAUSS_FILT_SZ_BY2 + 1] = {
  {0.046701777738928, 0.041214174199798, 0.028326060061745, 0.01516184737296,},
  {0.041214174199798, 0.036371381073904, 0.024997660266915, 0.013380283344101},
  {0.028326060061745, 0.024997660266915, 0.017180623896310, 0.009196125289586},
  {0.015161847372964, 0.013380283344101, 0.009196125289586, 0.004922331159344},
};

int HarrisCornerDetect(unsigned char *buf, int width, int height, int stride,
                       int *points, int max_points) {
  const double eps = 1e-6;
  int i, j, num_points = 0;
  double *gx2 = (double *)malloc(width * height * sizeof(double));
  double *gy2 = (double *)malloc(width * height * sizeof(double));
  double *gxy = (double *)malloc(width * height * sizeof(double));
  double *met = (double *)malloc(width * height * sizeof(double));
  double max_met = 0;
  double thresh;
  if (width < GAUSS_FILT_SZ + 2 || height < GAUSS_FILT_SZ + 2)
    return 0;
  // Sobel operator
  for (i = 1; i < height - 1; ++i) {
    for (j = 1; j < width - 1; ++j) {
      gx2[i * width + j] = (buf[(i + 1) * stride + j] - buf[(i - 1) * stride + j]) * 2 +
                           (buf[(i + 1) * stride + j - 1] - buf[(i - 1) * stride + j - 1]) +
                           (buf[(i + 1) * stride + j + 1] - buf[(i - 1) * stride + j + 1]);
      gy2[i * width + j] = (buf[i * stride + j + 1] - buf[i * stride + j - 1]) * 2 +
                           (buf[(i - 1) * stride + j + 1] - buf[(i - 1) * stride + j - 1]) +
                           (buf[(i + 1) * stride + j + 1] - buf[(i + 1) * stride + j - 1]);
      gx2[i * width + j] /= 8;
      gy2[i * width + j] /= 8;
      gxy[i * width + j] = gx2[i * width + j] * gy2[i * width + j];
      gx2[i * width + j] *= gx2[i * width + j];
      gy2[i * width + j] *= gy2[i * width + j];
    }
  }
  // Gaussian filter and metric computation
  for (i = 1 + GAUSS_FILT_SZ_BY2; i < height - 1 - GAUSS_FILT_SZ_BY2; ++i) {
    for (j = 1 + GAUSS_FILT_SZ_BY2; j < width - 1 - GAUSS_FILT_SZ_BY2; ++j) {
      double valx2 = 0;
      double valy2 = 0;
      double valxy = 0;
      int k, l;
      valx2 = GaussianKernel[0][0] * gx2[i * width + j];
      valy2 = GaussianKernel[0][0] * gy2[i * width + j];
      valxy = GaussianKernel[0][0] * gxy[i * width + j];
      for (k = 1; k <= GAUSS_FILT_SZ_BY2; k++) {
        valx2 += GaussianKernel[k][0] * (gx2[(i + k) * width + j] +
                                         gx2[(i - k) * width + j]);
        valy2 += GaussianKernel[k][0] * (gy2[(i + k) * width + j] +
                                         gy2[(i - k) * width + j]);
        valxy += GaussianKernel[k][0] * (gxy[(i + k) * width + j] +
                                         gxy[(i - k) * width + j]);
      }
      for (l = 1; l <= GAUSS_FILT_SZ_BY2; l++) {
        valx2 += GaussianKernel[0][l] * (gx2[i * width + j + l] +
                                         gx2[i * width + j - l]);
        valy2 += GaussianKernel[0][l] * (gy2[i * width + j + l] +
                                         gy2[i * width + j - l]);
        valxy += GaussianKernel[0][l] * (gxy[i * width + j + l] +
                                         gxy[i * width + j - l]);
      }
      for (k = 1; k <= GAUSS_FILT_SZ_BY2; k++) {
        for (l = 1; l <= GAUSS_FILT_SZ_BY2; ++l) {
          valx2 += GaussianKernel[k][l] * (gx2[(i + k) * width + (j + l)] +
                                           gx2[(i - k) * width + (j + l)] +
                                           gx2[(i + k) * width + (j - l)] +
                                           gx2[(i - k) * width + (j - l)]);
          valy2 += GaussianKernel[k][l] * (gy2[(i + k) * width + (j + l)] +
                                           gy2[(i - k) * width + (j + l)] +
                                           gy2[(i + k) * width + (j - l)] +
                                           gy2[(i - k) * width + (j - l)]);
          valxy += GaussianKernel[k][l] * (gxy[(i + k) * width + (j + l)] +
                                           gxy[(i - k) * width + (j + l)] +
                                           gxy[(i + k) * width + (j - l)] +
                                           gxy[(i - k) * width + (j - l)]);
        }
      }
      met[i * width + j] = (valx2 * valy2 - valxy * valxy)/(valx2 + valy2 + eps);
      if (met[i * width + j] > max_met) max_met = met[i * width + j];
    }
  }
  thresh = max_met < 2*HARRIS_THRESH ? max_met/2 : HARRIS_THRESH;
  free(gx2);
  free(gy2);
  free(gxy);
  for (i = 1 + GAUSS_FILT_SZ_BY2; i < height - 1 - GAUSS_FILT_SZ_BY2; ++i) {
    for (j = 1 + GAUSS_FILT_SZ_BY2; j < width - 1 - GAUSS_FILT_SZ_BY2; ++j) {
      if (met[i * width + j] > thresh) {
        points[num_points * 2] = j;
        points[num_points * 2 + 1] = i;
        num_points++;
        if (num_points == max_points) {
          free(met);
          return num_points;
        }
      }
    }
  }
  free(met);
  return NonmaxSuppression(buf, width, height, stride, points, num_points);
}
