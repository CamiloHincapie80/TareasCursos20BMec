"""Vista de la órbita baja, el impulso y la cónica heliocéntrica."""
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt6 import QtCore, QtWidgets

from .constantes import AU, DAY, RAD
from .hohmann import R1, R2, R3
from .transferencia_leo import (
    EARTH_RADIUS, EARTH_SOI, JUPITER_PHASE, JUPITER_RATE, circle, hohmann_burn,
    nominal_saturn_leg, nominal_scattering, transfer,
)


class TransferView(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        title = QtWidgets.QLabel('Transferencia Tierra–Júpiter · 20 de agosto de 1977, 14:29 UTC')
        title.setObjectName('pageTitle')
        layout.addWidget(title)
        controls = QtWidgets.QHBoxLayout()
        layout.addLayout(controls)
        self.altitude = QtWidgets.QDoubleSpinBox(suffix=' km')
        self.altitude.setRange(100, 2000)
        self.altitude.setValue(500)
        self.burn = QtWidgets.QDoubleSpinBox(suffix=' km/s')
        self.burn.setDecimals(8)
        self.burn.setSingleStep(.05)
        self.burn.setRange(0, 20)
        for label, widget in [('Altura de la órbita circular', self.altitude),
                              ('Impulso tangencial Δv₀', self.burn)]:
            controls.addWidget(QtWidgets.QLabel(label))
            controls.addWidget(widget)
            widget.setMinimumWidth(125)
        nominal = QtWidgets.QPushButton('Usar Δv de Hohmann')
        nominal.clicked.connect(self.reset_burn)
        controls.addWidget(nominal)
        controls.addStretch()
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.figure = Figure(figsize=(11, 6), layout='constrained', facecolor='#0b1020')
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas, 1)
        playback = QtWidgets.QHBoxLayout()
        layout.addLayout(playback)
        self.play = QtWidgets.QPushButton('Recorrer transferencia')
        self.play.clicked.connect(self.toggle_playback)
        playback.addWidget(self.play)
        self.progress = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.progress.setRange(0, 1000)
        self.progress.valueChanged.connect(self.place_marker)
        playback.addWidget(self.progress, 1)
        self.stage = QtWidgets.QLabel()
        self.stage.setMinimumWidth(155)
        playback.addWidget(self.stage)
        self.timer = QtCore.QTimer(self, interval=45)
        self.timer.timeout.connect(self.advance)
        note = QtWidgets.QLabel(
            'La vista solar usa órbitas circulares, con Júpiter en la fase de Hohmann. '
            'El empalme usa el exceso asintótico de la salida terrestre. La fecha fija el '
            'instante del impulso. En Júpiter se empalman las velocidades asintóticas y Saturno tiene la fase del diseño nominal. El control recorre la geometría por etapas.'
        )
        note.setWordWrap(True)
        note.setObjectName('muted')
        layout.addWidget(note)
        self.altitude.valueChanged.connect(self.redraw)
        self.burn.valueChanged.connect(self.redraw)
        self.reset_burn()

    def toggle_playback(self):
        if self.timer.isActive():
            self.stop_playback()
            return
        if self.progress.value() == self.progress.maximum():
            self.progress.setValue(0)
        self.play.setText('Pausar')
        self.timer.start()

    def stop_playback(self):
        self.timer.stop()
        self.play.setText('Recorrer transferencia')

    def advance(self):
        self.progress.setValue(self.progress.value()+8)
        if self.progress.value() == self.progress.maximum():
            self.stop_playback()

    def place_marker(self):
        remaining = self.progress.value()/self.progress.maximum()*len(self.stages)
        current = min(int(remaining),len(self.stages)-1)
        for i,(name,path,marker) in enumerate(self.stages):
            if i == current:
                point = path[round((remaining-current)*(len(path)-1))]
                marker.set_data([point[0]],[point[1]])
                self.stage.setText(name)
            marker.set_visible(i==current)
        self.canvas.draw_idle()

    def reset_burn(self):
        self.burn.blockSignals(True)
        self.burn.setValue(hohmann_burn(self.altitude.value()))
        self.burn.blockSignals(False)
        self.redraw()

    def redraw(self):
        altitude = self.altitude.value()
        chosen = self.burn.value()
        result = transfer(altitude, chosen)
        self.stop_playback()
        status = f'Circular {result.circular_speed:.3f} km/s · Después del impulso {result.burn_speed:.3f} km/s'
        status += f' · Δv para Hohmann a esta altura: {hohmann_burn(altitude):.6f} km/s'
        if result.earth_excess is None:
            arrival = 'El impulso queda por debajo del escape terrestre de dos cuerpos.'
        elif result.arrival_days is None:
            arrival = f'Exceso terrestre {result.earth_excess:.3f} km/s. La cónica solar queda dentro de la órbita de Júpiter.'
        else:
            arrival = (f'Exceso terrestre {result.earth_excess:.3f} km/s. Cruce del radio de Júpiter '
                       f'a los {result.arrival_days:.2f} días · desfase angular {result.phase_error_deg:+.2f}°.')
            if result.nominal:
                arrival += ' Llegada tangente al afelio de Hohmann.'
        scattering = nominal_scattering() if result.nominal else None
        if scattering is not None:
            saturn_leg = nominal_saturn_leg()
            saturn_miss = np.linalg.norm(saturn_leg.periapsis_position-saturn_leg.saturn_at_periapsis)
            arrival += (f'\nDispersión joviana: rₚ = {scattering.periapsis_km:,.0f} km '
                        f'· δ = {np.degrees(scattering.turn_rad):.2f}° '
                        f'· |v∞| = {np.linalg.norm(scattering.incoming):.3f} km/s '
                        f'· Δε solar = {scattering.energy_change:+.2f} km²/s². '
                        f'Saturno: {saturn_miss:,.0f} km a los {saturn_leg.days_to_periapsis:.1f} días desde Júpiter; se muestran {saturn_leg.days_after_periapsis:.0f} días posteriores.')
        elif result.arrival_days is not None:
            arrival += ' El sobrevuelo requiere ajustar la geometría de llegada.'
        self.summary.setText(status+'\n'+arrival)

        self.figure.clear()
        grid = self.figure.add_gridspec(2, 3 if scattering else 2,
                                        width_ratios=(1,1,2) if scattering else (1,2))
        leo = self.figure.add_subplot(grid[0, 0])
        escape = self.figure.add_subplot(grid[1, 0])
        solar = self.figure.add_subplot(grid[:, -1])
        plots = [(leo,'Órbita circular e impulso','km'), (escape,'Salida en el marco de la Tierra','km'),
                 (solar,'Transferencia en el marco solar','UA')]
        if scattering is not None:
            jovian = self.figure.add_subplot(grid[0, 1])
            saturn_view = self.figure.add_subplot(grid[1, 1])
            plots.extend([(jovian,'Dispersión en el marco de Júpiter','Millones de km'),
                          (saturn_view,'Sobrevuelo en el marco de Saturno','Millones de km')])
        for ax, title, units in plots:
            ax.set_title(title, color='#e6edf8')
            ax.set_facecolor('#121c30')
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel(units, color='#9eafc9')
            ax.tick_params(colors='#9eafc9', labelsize=8)
            ax.grid(alpha=.14)

        angles = np.linspace(0, 2*np.pi, 361)
        leo.fill(*circle(EARTH_RADIUS, angles).T, color='#1f507a', label='Tierra')
        leo.plot(*result.leo.T, color='#55b7fc', lw=1.8, label='Órbita circular')
        radius = EARTH_RADIUS+altitude
        leo.set_xlim(-1.25*radius, 1.55*radius)
        leo.set_ylim(-1.2*radius, 1.2*radius)
        leo.scatter(radius, 0, color='#ffd15a', zorder=4)
        if chosen > 0:
            leo.annotate('', xy=(radius,.55*radius), xytext=(radius,0),
                         arrowprops={'arrowstyle':'->','color':'#ffd15a','lw':2})
        leo.annotate(f'Δv₀ = {chosen:.6f} km/s', (radius,0), xytext=(-5,-15),
                     textcoords='offset points', ha='right', color='#ffd15a', fontsize=8)
        escape.plot(*circle(EARTH_SOI, angles).T, '--', color='#677c9c', label='Esfera terrestre')
        escape.plot(*result.escape.T, color='#48d8ff', lw=1.8, label='Cónica tras el impulso')
        escape.scatter(0, 0, s=16, color='#55b7fc', zorder=4)
        if np.max(np.linalg.norm(result.escape, axis=1)) < EARTH_SOI/20:
            limit = 1.25*np.max(np.linalg.norm(result.escape, axis=1))
            escape.set_xlim(-limit,limit)
            escape.set_ylim(-limit,limit)

        for radius, color, name in [(R1,'#55b7fc','Órbita de la Tierra'),(R2,'#edb875','Órbita de Júpiter')]:
            solar.plot(*circle(radius/AU, angles).T, color=color, alpha=.55, lw=1.1, label=name)
        solar.scatter(0, 0, s=38, color='#ffd15a')
        solar.scatter(R1/AU, 0, color='#55b7fc', zorder=4)
        if result.solar is not None:
            path = result.solar/AU
            solar.plot(*path.T, color='#48d8ff', lw=2.2,
                       label='Elipse de Hohmann' if result.nominal else 'Trayectoria con el Δv elegido')
            if result.arrival_days is not None:
                phase = JUPITER_PHASE+JUPITER_RATE*result.arrival_days*DAY
                solar.scatter(*circle(R2/AU, [phase])[0], color='#edb875', zorder=5,
                              label='Júpiter al cruce')
                solar.scatter(*path[-1], color='#48d8ff', s=14, zorder=6)
        if scattering is not None:
            hyperbola = scattering.hyperbola/1e6
            jovian.fill(*circle(RAD['J']/1e6, angles).T, color='#edb875')
            jovian.plot(*hyperbola.T, color='#48d8ff', lw=2)
            mid = len(hyperbola)//2
            jovian.scatter(*hyperbola[mid], s=16, color='#ffd15a', zorder=4)
            for vector, anchor, color in [(scattering.incoming,hyperbola[175],'#48d8ff'),
                                          (scattering.outgoing,hyperbola[-176],'#8ee0a4')]:
                direction = vector/np.linalg.norm(vector)*1.25
                jovian.annotate('', xy=anchor+direction, xytext=anchor,
                                arrowprops={'arrowstyle':'->','color':color,'lw':1.5})
            solar.plot(*circle(R3/AU, angles).T, color='#d4ca8c', alpha=.55, lw=1.1, label='Órbita de Saturno')
            solar.plot(*(saturn_leg.approach/AU).T, color='#8ee0a4', lw=2, label='Salida hacia Saturno')
            solar.plot(*(saturn_leg.departure/AU).T, color='#feb778', lw=2, label='Después de Saturno')
            solar.scatter(*(saturn_leg.saturn_at_periapsis/AU), color='#d4ca8c', zorder=5)
            saturn_track = saturn_leg.relative_path/1e6
            saturn_mid = len(saturn_track)//2
            for segment, color in [(saturn_leg.encounter_path[:saturn_mid+1],'#8ee0a4'),
                                   (saturn_leg.encounter_path[saturn_mid:],'#feb778')]:
                solar.plot(*(segment/AU).T, color=color, lw=2)
            saturn_view.fill(*circle(RAD['S']/1e6, angles).T, color='#d4ca8c')
            saturn_view.plot(*saturn_track[:saturn_mid+1].T, color='#8ee0a4', lw=2)
            saturn_view.plot(*saturn_track[saturn_mid:].T, color='#feb778', lw=2)
            saturn_view.scatter(*saturn_track[saturn_mid], color='#ffd15a', s=15, zorder=5)
        for ax in (escape, solar):
            ax.legend(loc='upper left', facecolor='#121c30', edgecolor='#25334c',
                      labelcolor='#e6edf8', fontsize=8)
        stages = [('Órbita circular', result.leo, leo),
                  ('Cónica tras el impulso', result.escape, escape)]
        if result.solar is not None:
            stages.append(('Transferencia solar', result.solar/AU, solar))
        if scattering is not None:
            stages.extend([('Dispersión joviana', hyperbola, jovian),
                           ('Salida hacia Saturno', saturn_leg.approach/AU, solar),
                           ('Sobrevuelo de Saturno', saturn_track, saturn_view),
                           ('Después de Saturno', saturn_leg.departure/AU, solar)])
        self.stages = [(name, points, ax.plot([],[],marker='o',ms=6,color='white',zorder=7)[0])
                       for name,points,ax in stages]
        self.progress.setValue(0)
        self.place_marker()
