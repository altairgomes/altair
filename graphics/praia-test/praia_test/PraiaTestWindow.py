# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
### BEGIN LICENSE
# This file is in the public domain
### END LICENSE

from locale import gettext as _

from gi.repository import Gtk # pylint: disable=E0611
import logging
logger = logging.getLogger('praia_test')

from praia_test_lib import Window
from praia_test.AboutPraiaTestDialog import AboutPraiaTestDialog
from praia_test.PreferencesPraiaTestDialog import PreferencesPraiaTestDialog
from praia_test.HextractionDialog import HextractionDialog

# See praia_test_lib.Window.py for more details about how this class works
class PraiaTestWindow(Window):
    __gtype_name__ = "PraiaTestWindow"
    
    def finish_initializing(self, builder): # pylint: disable=E1002
        """Set up the main window"""
        super(PraiaTestWindow, self).finish_initializing(builder)

        self.AboutDialog = AboutPraiaTestDialog
        self.PreferencesDialog = PreferencesPraiaTestDialog

        # Code for other initialization actions should be added here.
    def on_buttonok_clicked(self, widget):
        title = self.ui.select1.get_active_text()
        if title == 'Header Extraction':
            Hextract = HextractionDialog()
            Hextract.run()
        else:
            print title
