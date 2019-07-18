#!/usr/bin/env python
# coding: utf-8
import numpy as np
import numpy.linalg as LA
from Util import MSE
from MotionEST import MotionEST
"""Exhaust Search:"""


class Exhaust(MotionEST):
  """
    Constructor:
        cur_f: current frame
        ref_f: reference frame
        blk_sz: block size
        wnd_size: search window size
        metric: metric to compare the blocks distrotion
    """

  def __init__(self, cur_f, ref_f, blk_size, wnd_size, metric=MSE):
    self.name = 'exhaust'
    self.wnd_sz = wnd_size
    self.metric = metric
    super(Exhaust, self).__init__(cur_f, ref_f, blk_size)

  """
    search method:
        cur_r: start row
        cur_c: start column
    """

  def search(self, cur_r, cur_c):
    min_loss = self.dist(cur_r, cur_c, [0, 0], self.metric)
    cur_x = cur_c * self.blk_sz
    cur_y = cur_r * self.blk_sz
    ref_x = cur_x
    ref_y = cur_y
    #search all validate positions and select the one with minimum distortion
    for y in xrange(cur_y - self.wnd_sz, cur_y + self.wnd_sz):
      for x in xrange(cur_x - self.wnd_sz, cur_x + self.wnd_sz):
        if 0 <= x < self.width - self.blk_sz and 0 <= y < self.height - self.blk_sz:
          loss = self.dist(cur_r, cur_c, [y - cur_y, x - cur_x], self.metric)
          if loss < min_loss:
            min_loss = loss
            ref_x = x
            ref_y = y
    return ref_x, ref_y

  def est(self):
    for i in xrange(self.num_row):
      for j in xrange(self.num_col):
        ref_x, ref_y = self.search(i, j)
        self.mf[i, j] = np.array(
            [ref_y - i * self.blk_sz, ref_x - j * self.blk_sz])


"""Exhaust with Neighbor Constraint"""


class ExhaustNeighbor(MotionEST):
  """
    Constructor:
        cur_f: current frame
        ref_f: reference frame
        blk_sz: block size
        wnd_size: search window size
        beta: neigbor loss weight
        metric: metric to compare the blocks distrotion
    """

  def __init__(self, cur_f, ref_f, blk_size, wnd_size, beta, metric=MSE):
    self.name = 'exhaust + neighbor'
    self.wnd_sz = wnd_size
    self.beta = beta
    self.metric = metric
    super(ExhaustNeighbor, self).__init__(cur_f, ref_f, blk_size)
    self.assign = np.zeros((self.num_row, self.num_col), dtype=np.bool)

  """
    estimate neighbor loss:
        cur_r: current row
        cur_c: current column
        mv: current motion vector
    """

  def neighborLoss(self, cur_r, cur_c, mv):
    loss = 0
    #accumulate difference between current block's motion vector with neighbors'
    for i, j in {(-1, 0), (1, 0), (0, 1), (0, -1)}:
      nb_r = cur_r + i
      nb_c = cur_c + j
      if 0 <= nb_r < self.num_row and 0 <= nb_c < self.num_col and self.assign[
          nb_r, nb_c]:
        loss += LA.norm(mv - self.mf[nb_r, nb_c])
    return loss

  """
    search method:
        cur_r: start row
        cur_c: start column
    """

  def search(self, cur_r, cur_c):
    dist_loss = self.dist(cur_r, cur_c, [0, 0], self.metric)
    nb_loss = self.neighborLoss(cur_r, cur_c, np.array([0, 0]))
    min_loss = dist_loss + self.beta * nb_loss
    cur_x = cur_c * self.blk_sz
    cur_y = cur_r * self.blk_sz
    ref_x = cur_x
    ref_y = cur_y
    #search all validate positions and select the one with minimum distortion
    # as well as weighted neighbor loss
    for y in xrange(cur_y - self.wnd_sz, cur_y + self.wnd_sz):
      for x in xrange(cur_x - self.wnd_sz, cur_x + self.wnd_sz):
        if 0 <= x < self.width - self.blk_sz and 0 <= y < self.height - self.blk_sz:
          dist_loss = self.dist(cur_r, cur_c, [y - cur_y, x - cur_x],
                                self.metric)
          nb_loss = self.neighborLoss(cur_r, cur_c, [y - cur_y, x - cur_x])
          loss = dist_loss + self.beta * nb_loss
          if loss < min_loss:
            min_loss = loss
            ref_x = x
            ref_y = y
    return ref_x, ref_y

  def est(self):
    for i in xrange(self.num_row):
      for j in xrange(self.num_col):
        ref_x, ref_y = self.search(i, j)
        self.mf[i, j] = np.array(
            [ref_y - i * self.blk_sz, ref_x - j * self.blk_sz])
        self.assign[i, j] = True
