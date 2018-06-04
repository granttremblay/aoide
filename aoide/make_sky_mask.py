#!/usr/bin/env python
'''
Aoide | Reduction & Analysis of MUSE observations
-------------------------------------------------
Dr. Grant R. Tremblay | Harvard-Smithsonian Center for Astrophysics
grant.tremblay @ cfa.harvard.edu

See the README associated with this repository for documentation & examples.
'''


from __future__ import print_function

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits as pyfits

"""
Interactive Masking of a FITS-file. The FITS-file must be provided upon
creating a new instance. If no mask is provided, the routine will create one
from scratch. Otherwise, the supplied mask will be modified.
The masking porcess is carried out using the mouse: An area to mask is selected
by moving the mouse over it while pressing the left button. To unmask an area,
use the right button. The cuts might be changed by clicking on the wheel.
Note that plt.show() must be called after an instance of MaskFrame has been
created!

"""


class MaskFrame:
    """
    Initiate an instance
    """

    def __init__(self, image, mask_name, cuts=(0, 10), extension=0):
        fits_ima = pyfits.open(image)
        self.true_arr = fits_ima[extension].data
        if len(self.true_arr.shape) == 3:
            self.true_arr = self.true_arr[0, :]
        fits_ima.close()
        self.mask_name = mask_name
        self.extension = extension

        if os.path.exists(mask_name):
            self.in_mask = pyfits.open(mask_name, mode='update')
            self.mask = self.in_mask[0].data
        else:
            self.in_mask = None
            self.mask = np.zeros(self.true_arr.shape, dtype='Int16')

        self.plot_arr = self.true_arr + (self.mask * 1e9)

        self.lo_cut = cuts[0]
        self.hi_cut = cuts[1]

        self.fig = plt.figure(figsize=(8,8))
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title('LEFT: Mask | RIGHT: Unmask | Wheel: Change cuts')
        self.im = self.ax.imshow(
            self.true_arr, origin='lower', interpolation='nearest', cmap='magma')

        self.update()

        self.xM = []
        self.yM = []

        self._connect()

    """
    Connect the button_***_events to the corresponding methods
    """

    def _connect(self):
        self.ax.figure.canvas.mpl_connect('button_press_event', self.__on)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.__off)

    """
    The actions that are carried out when a mouse button is pressed:
    """

    def __on(self, event):
        if event.button == 2:
            print('Current cut levels are: {}, {}'.format(
                self.lo_cut, self.hi_cut))
            new_c = input('Enter new cut levels as low,high e.g. 0,20: ')
            self.lo_cut = float(new_c.split(',')[0])
            self.hi_cut = float(new_c.split(',')[1])
            self.update()
        else:
            if event.inaxes != self.ax.axes:
                print('Out of bounds!')
                return
            self.xM.append(int(round(event.xdata)))
            self.yM.append(int(round(event.ydata)))

    """
    The actions that are carried out when a mouse button is released.
    """

    def __off(self, event):
        if event.inaxes != self.ax.axes:
            print('Out of bounds!')
            return
        else:
            self.xM.append(int(round(event.xdata)))
            self.yM.append(int(round(event.ydata)))

            if len(self.xM) == 2:
                if event.button == 1:
                    self.mask[min(self.yM):max(self.yM) + 1,
                              min(self.xM):max(self.xM) + 1] = 1
                elif event.button == 3:
                    self.mask[min(self.yM):max(self.yM) + 1,
                              min(self.xM):max(self.xM) + 1] = 0

                self.plot_arr = self.true_arr + (self.mask * 1e9)
                self.update()

        self.xM = []
        self.yM = []

    """
    This method updates the graphical interface:
    """

    def update(self):
        self.im.set_data(self.plot_arr[:, :])
        self.im.set_clim(vmin=self.lo_cut, vmax=self.hi_cut)
        self.im.axes.figure.canvas.draw()

    """
    Save the mask under the filename specified in FrameMask.__init__
    Note that unlike the other methods, this method must be called explicitely
    """

    def save_mask(self):
        extension = self.extension
        if self.in_mask == None:
            maskHDU = pyfits.PrimaryHDU(self.mask)
            maskHDU.writeto(self.mask_name, overwrite=True)
        else:
            self.in_mask[0].data = self.mask
            self.in_mask.flush()



def main():
    if len(sys.argv) == 3:
        make_mask = MaskFrame(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        make_mask = MaskFrame(
            sys.argv[1], sys.argv[2], extension=int(sys.argv[3]))
    plt.show()
    make_mask.save_mask()


if __name__ == '__main__':
    main()
