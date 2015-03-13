# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
### BEGIN LICENSE
# This file is in the public domain
### END LICENSE

from gi.repository import Gtk # pylint: disable=E0611

from praia_test_lib.helpers import get_builder

import gettext
from gettext import gettext as _
gettext.textdomain('praia-test')

class HextractionDialog(Gtk.Dialog):
    __gtype_name__ = "HextractionDialog"
    

    def __new__(cls):
        """Special static method that's automatically called by Python when 
        constructing a new instance of this class.
        
        Returns a fully instantiated HextractionDialog object.
        """
        builder = get_builder('HextractionDialog')
        new_object = builder.get_object('hextraction_dialog')
        new_object.finish_initializing(builder)
        return new_object

    def finish_initializing(self, builder):
        """Called when we're finished initializing.

        finish_initalizing should be called after parsing the ui definition
        and creating a HextractionDialog object with it in order to
        finish initializing the start of the new HextractionDialog
        instance.
        """
        # Get a reference to the builder and set up the signals.
        self.builder = builder
        self.ui = builder.get_ui(self)
        a = ['0', '1']
        b = ['Teste1', 'Teste2']
        for idx, value in enumerate(a):
            self.ui.head1.append(a[idx], b[idx])
        self.ui.head1.set_active(1)

    def on_btn_ok_clicked(self, widget, data=None):
        """The user has elected to save the changes.

        Called before the dialog returns Gtk.ResponseType.OK from run().
        """
        arqoutput = self.ui.entry1.get_text()
        header = self.ui.head1.get_active_text()
        offset = self.ui.entry2.get_text()
        diretorio = self.ui.pasta1.get_filename()
        subdir = self.ui.subdir1.get_active()
        checkdat = self.ui.checkdat1.get_active()
        print 'Vai rodar PRAIA'
        

    def on_btn_cancel_clicked(self, widget, data=None):
        """The user has elected cancel changes.

        Called before the dialog returns Gtk.ResponseType.CANCEL for run()
        """
        self.destroy()

    def on_btn_add_hd_clicked(self, widget, data=None):
        """The user has elected cancel changes.

        Called before the dialog returns Gtk.ResponseType.CANCEL for run()
        """
        print 'Vai adicionar header'
        
    def on_readme_clicked(self, widget, data=None):
        """The user has elected cancel changes.

        Called before the dialog returns Gtk.ResponseType.CANCEL for run()
        """
        print 'Vai mostrar README'


if __name__ == "__main__":
    dialog = HextractionDialog()
    dialog.show()
    Gtk.main()
