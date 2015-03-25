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

from gi.repository import Gtk # pylint: disable=E0611

from altair_lib.helpers import get_builder
from altair_lib.observation import Observation
from astropy.time import Time, TimeDelta

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
        
        self.def_tooltips()
        
        self.ui.horain.set_range(12,23)
        self.ui.horain.set_value(18)
        self.ui.horain.set_increments(1, 1)
        self.ui.horainmin.set_range(0,59)
        self.ui.horainmin.set_value(0)
        self.ui.horainmin.set_increments(1, 1)
        
        self.ui.horafin.set_range(0,12)
        self.ui.horafin.set_value(6)
        self.ui.horafin.set_increments(1, 1)
        self.ui.horafinmin.set_range(0,59)
        self.ui.horafinmin.set_value(0)
        self.ui.horafinmin.set_increments(1, 1)
        
        self.ui.intdata.set_range(0,120)
        self.ui.intdata.set_value(60)
        self.ui.intdata.set_increments(5, 20)

    def def_tooltips(self):
        self.ui.inputfile.set_tooltip_markup('Arquivo com os dados de entrada')
        self.ui.col_name.set_tooltip_markup('Números das colunas que serao os nomes dos objetos separados por virgula')
        self.ui.col_coord.set_tooltip_markup('Números das colunas que contem as coordenadas dos objetos separados por virgula')
        self.ui.col_comment.set_tooltip_markup('Números das colunas que serao comentarios dos objetos separados por virgula')
        self.ui.fuso.set_tooltip_markup('Fuso local')
        self.ui.latitude.set_tooltip_markup('Latitude do local')
        self.ui.longitude.set_tooltip_markup('Longitude do local')
        self.ui.altitude.set_tooltip_markup('Altitude do local')
        self.ui.limalt.set_tooltip_markup('Altura minima de observacao')
        self.ui.limdist.set_tooltip_markup('Tamanho do campo: Objetos proximos serao identificados')
        self.ui.diainicio.set_tooltip_markup('Dia do inicio da observacao')
        self.ui.horain.set_tooltip_markup('Hora do inicio da observacao')
        self.ui.horainmin.set_tooltip_markup('Minuto do inicio da observacao')
        self.ui.diafim.set_tooltip_markup('Dia do fim da observacao')
        self.ui.horafin.set_tooltip_markup('Hora do fim da observacao')
        self.ui.horafinmin.set_tooltip_markup('Minuto do fim da observacao')
        self.ui.intdata.set_tooltip_markup('Passo em que os dados serao mostrados')
        

    def on_btn_ok_clicked(self, widget, data=None):
        """The user has elected to save the changes.

        Called before the dialog returns Gtk.ResponseType.OK from run().
        """
        arqin = self.ui.inputfile.get_filename()
        col_name = self.ui.col_name.get_text()
        col_coord = self.ui.col_coord.get_text()
        col_comment = self.ui.col_comment.get_text()
        fuso = int(self.ui.fuso.get_text())
        latitude = self.ui.latitude.get_text()
        longitude = self.ui.longitude.get_text()
        altitude = int(self.ui.altitude.get_text())
        limalt = int(self.ui.limalt.get_text())
        limdist = int(self.ui.limdist.get_text())
        anoin, mesin, diain = self.ui.diainicio.get_date()
        horain = self.ui.horain.get_value_as_int()
        horainmin = self.ui.horainmin.get_value_as_int()
        anofin, mesfin, diafin = self.ui.diafim.get_date()
        horafin = self.ui.horafin.get_value_as_int()
        horafinmin = self.ui.horafinmin.get_value_as_int()
        intval = int(self.ui.intdata.get_value_as_int())
        
        cols_n, cols_c, cols_cm = None, None, None
        if len(col_name) != 0:
            cols_n = map(int, col_name.split(','))
        if len(col_coord) != 0:
            cols_c = map(int, col_coord.split(','))
        if len(col_comment) != 0:
            cols_cm = map(int, col_comment.split(','))
        
        obs = Observation(fuso, latitude, longitude, altitude)
        objs, coords, comments = obs.read(arqin, cols_n, cols_c, cols_cm)
        obs.create_plan(coords,objs, comments, '{}-{}-{} {}:{}:00.000'.format(anoin,mesin,diain,horain,horainmin), '{}-{}-{} {}:{}:00.000'.format(anofin,mesfin,diafin,horafin,horafinmin),
intval, limalt=limalt, size=limdist, path='/home/altairjunior/Documentos', uipg=self.ui.progressbar)

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
