# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
### BEGIN LICENSE
# Copyright (C) 2015 Altair Ramos altaigomesjr@gmail.com
# This program is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License version 3, as published 
# by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranties of 
# MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR 
# PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with this program.  If not, see <http://www.gnu.org/licenses/>.
### END LICENSE

from locale import gettext as _

from gi.repository import Gtk # pylint: disable=E0611
import logging
logger = logging.getLogger('altair')

from altair_lib import Window
from altair.AboutAltairDialog import AboutAltairDialog
from altair.PreferencesAltairDialog import PreferencesAltairDialog
from altair.MapaDialog import MapaDialog
from altair.ObservationDialog import ObservationDialog

# See altair_lib.Window.py for more details about how this class works
class AltairWindow(Window):
    __gtype_name__ = "AltairWindow"
    
    def finish_initializing(self, builder): # pylint: disable=E1002
        """Set up the main window"""
        super(AltairWindow, self).finish_initializing(builder)

        self.AboutDialog = AboutAltairDialog
        self.PreferencesDialog = PreferencesAltairDialog

        # Code for other initialization actions should be added here.
    def on_buttonok_clicked(self, widget):
        title = self.ui.select1.get_active_text()
        if title == 'Mapa':
            Mapa = MapaDialog()
            result = Mapa.run()
        elif title == 'Observation':
            Obs = ObservationDialog()
            result = Obs.run()
        else:
            print title
