# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
### BEGIN LICENSE
# This file is in the public domain
### END LICENSE

from locale import gettext as _

from gi.repository import Gtk # pylint: disable=E0611
import logging
logger = logging.getLogger('altair')

from altair_lib import Window
from altair.AboutAltairDialog import AboutAltairDialog
from altair.PreferencesAltairDialog import PreferencesAltairDialog
from altair.MapaDialog import MapaDialog

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
        else:
            print title
