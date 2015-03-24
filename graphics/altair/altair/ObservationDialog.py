# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
### BEGIN LICENSE
# This file is in the public domain
### END LICENSE

from gi.repository import Gtk # pylint: disable=E0611

from altair_lib.helpers import get_builder
from altair_lib.observation import Observation

import gettext
from gettext import gettext as _
gettext.textdomain('altair')

class ObservationDialog(Gtk.Dialog):
    __gtype_name__ = "ObservationDialog"

    def __new__(cls):
        """Special static method that's automatically called by Python when 
        constructing a new instance of this class.
        
        Returns a fully instantiated ObservationDialog object.
        """
        builder = get_builder('ObservationDialog')
        new_object = builder.get_object('observation_dialog')
        new_object.finish_initializing(builder)
        return new_object

    def finish_initializing(self, builder):
        """Called when we're finished initializing.

        finish_initalizing should be called after parsing the ui definition
        and creating a ObservationDialog object with it in order to
        finish initializing the start of the new ObservationDialog
        instance.
        """
        # Get a reference to the builder and set up the signals.
        self.builder = builder
        self.ui = builder.get_ui(self)
        self.ui.horain.set_range(12.0,24.0)
        self.ui.horain.set_value(18.0)
        self.ui.horain.set_increments(0.5, 1.0)
        
        self.ui.horafin.set_range(0.0,12.0)
        self.ui.horafin.set_value(6.0)
        self.ui.horafin.set_increments(0.5, 1.0)
        
        self.ui.intdata.set_range(0,120)
        self.ui.intdata.set_value(60)
        self.ui.intdata.set_increments(5, 20)

        

    def on_btn_ok_clicked(self, widget, data=None):
        """The user has elected to save the changes.

        Called before the dialog returns Gtk.ResponseType.OK from run().
        """
        arqin = self.ui.inputfile.get_filename()
        col_name = self.ui.col_name.get_text()
        col_coord = self.ui.col_coord.get_text()
        col_comment = self.ui.col_comment.get_text()
        fuso = self.ui.fuso.get_text()
        latitude = self.ui.latitude.get_text()
        longitude = self.ui.longitude.get_text()
        altitude = self.ui.altitude.get_text()
        limalt = self.ui.limalt.get_text()
        limdist = self.ui.limdist.get_text()
        
        cols_n = map(int, col_name.split(','))
        cols_c = map(int, col_coord.split(','))
        cols_cm = map(int, col_comment.split(','))
        obs = Observation(fuso, latitude, longitude, altitude)
        objs, coords, comments = observation.read(arqin, cols_n, cols_c, cols_cm)
        observation.create_plan(coords,objs, comments, horain, horafin, intinf, limalt=limalt, size=limdist)

    def on_btn_cancel_clicked(self, widget, data=None):
        """The user has elected cancel changes.

        Called before the dialog returns Gtk.ResponseType.CANCEL for run()
        """
        self.destroy()
        pass


if __name__ == "__main__":
    dialog = ObservationDialog()
    dialog.show()
    Gtk.main()
