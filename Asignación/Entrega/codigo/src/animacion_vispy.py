"""Interfaz local Qt con escenas VisPy 3D; datos calculados fuera del bucle gráfico."""
import json
import sys
import time
import numpy as np
from PyQt6 import QtCore,QtGui,QtWidgets
from vispy import app,scene
from .constantes import AU,DAY,MU
from .encuentros import calendar

FRAMES={'H':'Marco solar','T':'Marco de la Tierra','J':'Marco de Júpiter','S':'Marco de Saturno'}
TRAJECTORIES={'Hohmann':'Gravedad continua ajustada','K':'Esferas de influencia','Voyager':'Voyager 2'}

COLORS={'Voyager':'#b69bff','Hohmann':'#48d8ff','N':'#53cf9c','P':'#61a9ff','P0':'#d98cfa','K':'#ff6c7b',
        'T':'#55b7fc','J':'#edb875','S':'#d4ca8c','Sol':'#ffd15a'}
NAMES={'Voyager':'Voyager · Horizons original','Hohmann':'Hohmann + asistencia idealizada',
       'N':'N · referencia gravitacional','P':'P · cambio de marco exacto',
       'P0':'P0 · ecuación original','K':'K · cónicas empalmadas'}

def rgba(color,alpha=1.):
    return tuple(int(color[i:i+2],16)/255 for i in (1,3,5))+(alpha,)

def at(values,times,t):
    i=int(np.clip(np.searchsorted(times,t,side='right')-1,0,len(times)-1))
    if i==len(times)-1:return values[i],i
    frac=(t-times[i])/(times[i+1]-times[i])
    return values[i]*(1-frac)+values[i+1]*frac,i

class OrbitalPanel:
    def __init__(self,parent,times,positions,title,soi=None):
        self.times=times;self.positions=positions;self.soi=soi or {};self.origin='H'
        self.time=0.;self.show_soi=True;self.show_orbits=True;self.show_full=True
        self.canvas=scene.SceneCanvas(keys='interactive',bgcolor='#0b1220',size=(900,650),dpi=96)
        self.view=self.canvas.central_widget.add_view()
        self.view.camera=scene.cameras.TurntableCamera(fov=0,elevation=65,azimuth=-85,up='+z')
        self.title=title;self.full={};self.trail={};self.markers={};self.labels={};self.checked={}
        for key,values in positions.items():
            color=COLORS[key]
            self.full[key]=scene.visuals.Line(parent=self.view.scene,color=rgba(color,.55 if key in NAMES else .28),width=1.7 if key in NAMES else 1.2,method='gl',antialias=True)
            self.trail[key]=scene.visuals.Line(parent=self.view.scene,color=color,width=2.3 if key in NAMES else 1.2,method='gl',antialias=True)
            self.full[key].set_gl_state('translucent',depth_test=False)
            self.trail[key].set_gl_state('translucent',depth_test=False)
            self.markers[key]=scene.visuals.Markers(parent=self.view.scene)
            self.labels[key]=scene.visuals.Text({'T':'Tierra','J':'Júpiter','S':'Saturno'}.get(key,key),font_size=9,color=color,
                                                parent=self.view.scene,anchor_x='left',anchor_y='bottom')
            self.checked[key]=True
        self.sun=scene.visuals.Markers(parent=self.view.scene)
        self.sun_label=scene.visuals.Text('Sol',color=COLORS['Sol'],font_size=10,parent=self.view.scene,
                                         anchor_x='left',anchor_y='bottom')
        self.rings={k:[scene.visuals.Line(parent=self.view.scene,color=rgba(COLORS[k],.45),width=1.2,method='gl',antialias=True) for _ in range(3)] for k in self.soi}
        self.axis=scene.visuals.XYZAxis(parent=self.view.scene)
        from vispy.visuals.transforms import STTransform
        self.axis.transform=STTransform(scale=(.3,.3,.3))
        self.status=QtWidgets.QLabel(title);self.status.setWordWrap(True)
        container=QtWidgets.QWidget(parent);layout=QtWidgets.QVBoxLayout(container);layout.setContentsMargins(0,0,0,0)
        layout.addWidget(self.status)
        self.preview=None
        if isinstance(self.canvas.native,QtWidgets.QWidget):
            layout.addWidget(self.canvas.native,1)
        else:
            self.preview=QtWidgets.QLabel('Render EGL de prueba');self.preview.setMinimumSize(350,450)
            self.preview.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter);layout.addWidget(self.preview,1)
        self.widget=container
        self.update_geometry();self.reset_camera()

    def coords(self,key):
        # Restar en float64 antes de enviar float32 a OpenGL evita pérdida de precisión local.
        p=self.positions[key]
        if self.origin!='H':p=p-self.positions[self.origin]
        return np.asarray(p/AU,dtype=np.float32)

    def reset_camera(self,close=False):
        if self.origin=='H':
            arrays=[self.coords(key) for key in self.positions if self.checked[key]]
            p=np.concatenate(arrays) if arrays else np.zeros((1,3))
            p=p[np.all(np.isfinite(p),axis=1)]
            if not len(p):p=np.zeros((1,3))
            lo=p.min(axis=0);hi=p.max(axis=0)
            self.view.camera.center=tuple((hi+lo)/2)
            self.view.camera.scale_factor=max(2.,float(np.max(hi-lo))*1.25)
        else:
            self.view.camera.center=(0,0,0)
            radius=at(self.soi[self.origin],self.times,self.time)[0] if self.origin in self.soi else .3*AU
            self.view.camera.scale_factor=.03 if close else 3*float(radius)/AU
        self.view.camera.elevation=65;self.view.camera.azimuth=-85

    def update_geometry(self):
        for key in self.positions:
            p=self.coords(key)
            self.full[key].set_data(pos=p[np.all(np.isfinite(p),axis=1)])
        self.update(self.time)

    def set_origin(self,origin,close=False):
        self.origin=origin;self.update_geometry();self.reset_camera(close)

    def update(self,t):
        self.time=float(t)
        scale=float(self.view.camera.scale_factor or 10.)*.12
        self.axis.transform.scale=(scale,scale,scale)
        for key in self.positions:
            p=self.coords(key);point,i=at(p,self.times,t)
            available=np.all(np.isfinite(point))
            trail=p[:i+1]
            trail=trail[np.all(np.isfinite(trail),axis=1)]
            if available:
                trail=np.vstack((trail,point))
            self.trail[key].set_data(pos=trail)
            if not available:
                self.trail[key].visible=False;self.markers[key].visible=False;self.labels[key].visible=False
                self.full[key].visible=False
                continue
            self.markers[key].set_data(pos=np.asarray([point]),face_color=COLORS[key],edge_width=0,
                                       size=9 if key in NAMES else 13)
            self.labels[key].pos=point
            visible=self.checked[key]
            if key not in NAMES and not self.show_orbits:
                self.full[key].visible=False;self.trail[key].visible=False
            else:
                self.full[key].visible=visible and self.show_full
                self.trail[key].visible=visible
            self.markers[key].visible=visible;self.labels[key].visible=visible
        sun=np.zeros(3) if self.origin=='H' else -at(self.positions[self.origin],self.times,t)[0]/AU
        self.sun.set_data(pos=np.asarray([sun],dtype=np.float32),face_color=COLORS['Sol'],edge_width=0,size=18)
        self.sun_label.pos=sun
        ang=np.linspace(0,2*np.pi,181);zero=np.zeros_like(ang)
        for k,rings in self.rings.items():
            center=at(self.coords(k),self.times,t)[0]
            radius=float(at(self.soi[k],self.times,t)[0])/AU
            circles=[np.c_[np.cos(ang),np.sin(ang),zero],np.c_[np.cos(ang),zero,np.sin(ang)],np.c_[zero,np.cos(ang),np.sin(ang)]]
            for line,circle in zip(rings,circles):
                line.set_data(pos=np.asarray(center+radius*circle,dtype=np.float32));line.visible=self.show_soi and (self.origin=='H' or radius<self.view.camera.scale_factor)
        self.canvas.update()

STYLE='''
QWidget { background:#0b1020; color:#e6edf8; font-family:"DejaVu Sans"; font-size:12px; }
QFrame#card { background:#121c30; border:1px solid #25334c; border-radius:16px; }
QFrame#card QLabel, QFrame#card QCheckBox { background:transparent; }
QLabel#pageTitle { font-size:23px; font-weight:700; color:#f5f8ff; }
QLabel#cardTitle { font-size:19px; font-weight:700; }
QLabel#muted { color:#9eafc9; font-size:11px; }
QLabel#eyebrow { color:#78dcca; font-size:10px; font-weight:700; padding-top:6px; }
QLabel#date { color:#b9ccec; font-size:12px; }
QPushButton,QComboBox,QDoubleSpinBox { background:#1d2c45; border:1px solid #334766; border-radius:7px; padding:7px; }
QPushButton:hover { background:#2c4264; border-color:#70c8df; }
QPushButton#primary { background:#236b80; border-color:#398da2; font-weight:600; }
QPushButton#primary:hover { background:#2a8198; }
QPushButton[trajectory="true"] { color:#9eafc9; background:#121c30; font-weight:600; }
QPushButton[trajectory="true"]:checked { color:#f5f8ff; background:#2c4264; border-color:#70c8df; }
QSlider::groove:horizontal { height:5px; background:#293b56; border-radius:2px; }
QSlider::sub-page:horizontal { background:#67d5c5; border-radius:2px; }
QSlider::handle:horizontal { background:#c2fff2; width:13px; margin:-5px 0; border-radius:6px; }
QTableWidget { background:#0f192b; border:1px solid #263953; border-radius:8px; gridline-color:#223149; }
QHeaderView::section { background:#1b2940; color:#adbed8; border:none; padding:8px 2px; font-size:10px; }
QCheckBox { spacing:5px; }
QScrollArea { border:none; }
'''

class MetricCard(QtWidgets.QFrame):
    def __init__(self,key,summary,saturn,parent=None):
        super().__init__(parent);self.setObjectName('card');self.key=key
        layout=QtWidgets.QVBoxLayout(self);layout.setContentsMargins(18,14,18,14);layout.setSpacing(8)
        name=TRAJECTORIES[key]
        title=QtWidgets.QLabel('●  '+name);title.setObjectName('cardTitle');title.setStyleSheet('color:'+COLORS[key]);layout.addWidget(title)
        self.values={};grid=QtWidgets.QGridLayout();grid.setHorizontalSpacing(28);layout.addLayout(grid)
        for i,(field,label) in enumerate([('speed','Velocidad heliocéntrica'),('energy','Energía específica solar'),
                                         ('distance','Distancia mínima a Júpiter'),('angle','Ángulo de dispersión δ'),
                                         ('saturn_distance','Distancia mínima a Saturno'),('saturn_date','Encuentro con Saturno')]):
            caption=QtWidgets.QLabel(label);caption.setObjectName('muted')
            value=QtWidgets.QLabel();value.setStyleSheet('font-size:17px;font-weight:600;')
            cell=QtWidgets.QVBoxLayout();cell.addWidget(caption);cell.addWidget(value);grid.addLayout(cell,i//2,i%2)
            self.values[field]=value
        valid=not summary['boundary_minimum']
        self.values['distance'].setText(f"{summary['rp_km']:,.1f} km"+(' *' if not valid else ''))
        angle=summary.get('delta_deg')
        self.values['angle'].setText(f'{angle:.3f}°' if valid and angle is not None else '—')
        self.values['saturn_distance'].setText(f"{saturn['rp_km']:,.1f} km")
        self.values['saturn_date'].setText(calendar(saturn['jd_tdb'])[:10])
        self.values['saturn_date'].setToolTip(calendar(saturn['jd_tdb']))
        note=QtWidgets.QLabel(summary['description']);note.setWordWrap(True);note.setObjectName('muted');layout.addWidget(note)

    def update_state(self,position,velocity):
        if not np.all(np.isfinite(position)):
            self.values['speed'].setText('Sin datos');self.values['energy'].setText('Sin datos')
            return
        speed=float(np.linalg.norm(velocity));energy=.5*speed**2-MU/np.linalg.norm(position)
        self.values['speed'].setText(f'{speed:,.3f} km/s')
        self.values['energy'].setText(f'{energy:+,.3f} km²/s²')


class AnimationWindow(QtWidgets.QMainWindow):
    def __init__(self,data):
        super().__init__();self.data=data;self.meta=json.loads(str(data['metadata']))
        self.times=data['t'];self.duration=float(self.times[-1]);self.t=0.;self.playing=False;self.last_tick=time.monotonic()
        self.setWindowTitle('Orbital · Comparación de trayectorias');self.resize(1480,1000)
        central=QtWidgets.QWidget()
        root=QtWidgets.QVBoxLayout(central);root.setContentsMargins(22,18,22,18);root.setSpacing(10)
        from .vista_transferencia import TransferView
        self.tabs=QtWidgets.QTabWidget()
        self.transfer_view=TransferView(self.tabs)
        self.tabs.addTab(self.transfer_view,'Transferencia de Hohmann')
        self.tabs.addTab(central,'Comparación hasta Saturno')
        self.setCentralWidget(self.tabs)
        title=QtWidgets.QLabel('ORBITAL  /  Tres trayectorias, una misma escena');title.setObjectName('pageTitle');root.addWidget(title)
        subtitle=QtWidgets.QLabel('Marcos con ejes eclípticos J2000 · tiempo transcurrido compartido · distancias en UA · arrastra para girar y usa la rueda para acercar')
        subtitle.setObjectName('muted');root.addWidget(subtitle)
        legend=QtWidgets.QHBoxLayout();root.addLayout(legend)
        self.visibility={}
        for key,name in TRAJECTORIES.items():
            button=QtWidgets.QPushButton(name)
            button.setProperty('trajectory',True)
            button.setCheckable(True)
            button.setChecked(True)
            swatch=QtGui.QPixmap(10,10)
            swatch.fill(QtGui.QColor(COLORS[key]))
            button.setIcon(QtGui.QIcon(swatch))
            button.setToolTip('Mostrar u ocultar '+name)
            button.toggled.connect(lambda visible,k=key:self.set_visible(k,visible))
            legend.addWidget(button)
            self.visibility[key]=button
        legend.addStretch()
        self.frame=QtWidgets.QComboBox()
        for origin,label in FRAMES.items():self.frame.addItem(label,origin)
        legend.addWidget(self.frame)
        reset=QtWidgets.QPushButton('Cámara inicial');reset.clicked.connect(lambda:self.panel.reset_camera());legend.addWidget(reset)
        positions={k:data[k] for k in ['Hohmann','K','Voyager','T','J','S']}
        self.panel=OrbitalPanel(central,self.times,positions,'',soi={k:data['soi_'+k] for k in ('T','J','S')})
        self.panel.canvas.size=(1250,480);self.panel.status.setObjectName('muted')
        self.panel.labels['K'].text='Esferas';self.panel.labels['Hohmann'].text='Gravedad continua';self.panel.labels['Voyager'].text='Voyager 2'
        self.panel.widget.setMinimumHeight(270)
        if self.panel.preview:self.panel.preview.setMinimumSize(350,270)
        root.addWidget(self.panel.widget,1)
        self.frame.currentIndexChanged.connect(self.change_frame)
        self.panel.status.hide()
        self.date=QtWidgets.QLabel();self.date.setObjectName('date');root.addWidget(self.date)
        self.slider=QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal);self.slider.setRange(0,10000)
        self.slider.valueChanged.connect(lambda v:self.seek(self.duration*v/10000,False));root.addWidget(self.slider)
        controls=QtWidgets.QHBoxLayout();root.addLayout(controls)
        self.play=QtWidgets.QPushButton('▶  Reproducir');self.play.setObjectName('primary');self.play.clicked.connect(self.toggle);controls.addWidget(self.play)
        start=QtWidgets.QPushButton('↺  Inicio');start.clicked.connect(lambda:self.seek(0));controls.addWidget(start)
        controls.addWidget(QtWidgets.QLabel('Avance:'))
        self.speed=QtWidgets.QDoubleSpinBox();self.speed.setRange(.01,1000);self.speed.setValue(35);self.speed.setSuffix(' días/s');controls.addWidget(self.speed)
        self.loop=QtWidgets.QCheckBox('Repetir');controls.addWidget(self.loop)
        for label,attribute in [('Recorridos completos','show_full'),('Órbitas planetarias','show_orbits'),('Esferas de influencia','show_soi')]:
            box=QtWidgets.QCheckBox(label);box.setChecked(True)
            box.toggled.connect(lambda visible,a=attribute:self.option(a,visible));controls.addWidget(box)
        self.encounters=QtWidgets.QComboBox();self.encounters.addItem('Ir a un evento…',None)
        for key,name in [('Hohmann','Gravedad continua'),('K','Esferas'),('Voyager','Voyager')]:
            self.encounters.addItem('Saturno · '+name,('S',self.meta['saturn_encounters'][key]['elapsed_days']*DAY))
        for key,departure in self.meta['earth_departures'].items():
            if departure['status']=='crossing':
                self.encounters.addItem('Salida terrestre · '+TRAJECTORIES[key],('T',departure['elapsed_seconds']))
        self.encounters.addItem('Tierra · lanzamiento',('T',0.))
        self.encounters.currentIndexChanged.connect(self.jump);controls.addWidget(self.encounters)
        controls.addStretch()
        cards=QtWidgets.QHBoxLayout();cards.setSpacing(14);root.addLayout(cards)
        self.cards={}
        for key in ('Hohmann','K','Voyager'):
            card=MetricCard(key,self.meta['jupiter_encounters'][key],self.meta['saturn_encounters'][key],central);cards.addWidget(card,1);self.cards[key]=card
        foot=QtWidgets.QLabel(self.meta['comparison_note']);foot.setObjectName('muted');foot.setWordWrap(True);root.addWidget(foot)
        for key,callback in [('Space',self.toggle),('R',lambda:self.panel.reset_camera())]:
            shortcut=QtGui.QShortcut(QtGui.QKeySequence(key),central)
            shortcut.setContext(QtCore.Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(callback)
        self.tabs.currentChanged.connect(self.pause_hidden_view)
        self.timer=QtCore.QTimer(self);self.timer.timeout.connect(self.tick);self.timer.start(33);self.seek(0)

    def pause_hidden_view(self,index):
        if index==0 and self.playing:
            self.toggle()
        elif index==1:
            self.transfer_view.stop_playback()

    def change_frame(self,index):
        origin=self.frame.itemData(index)
        self.panel.set_origin(origin)
        self.panel.status.setVisible(origin=='T')
        if origin=='T':
            self.update_earth_status()

    def update_earth_status(self):
        radius=at(self.data['soi_T'],self.times,self.t)[0]
        earth=at(self.data['T'],self.times,self.t)[0]
        labels=[]
        for key,departure in self.meta['earth_departures'].items():
            distance=np.linalg.norm(at(self.data[key],self.times,self.t)[0]-earth)
            if not np.isfinite(distance):
                labels.append(f'{TRAJECTORIES[key]} sin datos en esta fecha')
                continue
            if departure['status']=='crossing':
                event='salida '+departure['date_tdb'].replace('T',' ',1)
            elif departure['status']=='no_data':
                event='cruce anterior a los datos'
            elif departure['status']=='outside_at_start':
                event='ya fuera al inicio'
            else:
                event='sin cruce en los primeros tres días'
            labels.append(f'{TRAJECTORIES[key]} {distance:,.0f} km ({event})')
        self.panel.status.setText(f'Esfera terrestre: {radius:,.0f} km · '+' | '.join(labels))

    def jump(self,index):
        target=self.encounters.itemData(index)
        if target is None:return
        if self.playing:self.toggle()
        origin,seconds=target
        self.seek(seconds)
        self.frame.setCurrentIndex(self.frame.findData(origin))
        self.panel.reset_camera(close=origin=='S')
        self.encounters.blockSignals(True);self.encounters.setCurrentIndex(0);self.encounters.blockSignals(False)

    def set_visible(self,key,visible):
        self.panel.checked[key]=visible;self.panel.update(self.t)

    def option(self,attribute,visible):
        setattr(self.panel,attribute,visible);self.panel.update(self.t)

    def seek(self,t,update_slider=True):
        self.t=float(np.clip(t,0,self.duration));self.last_tick=time.monotonic();self.panel.update(self.t)
        if self.panel.origin=='T':self.update_earth_status()
        if update_slider:
            self.slider.blockSignals(True);self.slider.setValue(round(10000*self.t/self.duration));self.slider.blockSignals(False)
        self.date.setText(f"Día {self.t/DAY:,.2f}  /  {self.duration/DAY:,.0f}   ·   "+calendar(self.meta['epoch_jd']+self.t/DAY))
        for key,card in self.cards.items():
            card.update_state(at(self.data[key],self.times,self.t)[0],at(self.data['v_'+key],self.times,self.t)[0])

    def toggle(self):
        if not self.playing and self.t>=self.duration:self.seek(0)
        self.playing=not self.playing;self.last_tick=time.monotonic();self.play.setText('Ⅱ  Pausar' if self.playing else '▶  Reproducir')

    def tick(self):
        now=time.monotonic();elapsed=min(now-self.last_tick,.25);self.last_tick=now
        if not self.playing:return
        t=self.t+elapsed*self.speed.value()*DAY
        if t>=self.duration:
            if self.loop.isChecked():t%=self.duration
            else:t=self.duration;self.toggle()
        self.seek(t)

    def closeEvent(self,event):
        self.timer.stop();self.panel.canvas.close();super().closeEvent(event)


def launch(data,smoke=False,headless=False):
    app.use_app('egl' if headless else 'pyqt6')
    qt=QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1]);qt.setStyle('Fusion');qt.setStyleSheet(STYLE)
    window=AnimationWindow(data);window.show()
    if smoke:
        def verify():
            try:
                view=window.transfer_view
                assert 'Llegada tangente' in view.summary.text()
                view.burn.setValue(1.)
                assert 'por debajo del escape' in view.summary.text()
                view.altitude.setValue(700.)
                view.reset_burn()
                assert 'Llegada tangente' in view.summary.text()
                view.progress.setValue(200)
                assert view.stage.text()=='Cónica tras el impulso'
                view.progress.setValue(500)
                assert view.stage.text()=='Dispersión joviana'
                view.progress.setValue(650)
                assert view.stage.text()=='Salida hacia Saturno'
                view.progress.setValue(800)
                assert view.stage.text()=='Sobrevuelo de Saturno'
                view.progress.setValue(950)
                assert view.stage.text()=='Después de Saturno'
                window.tabs.setCurrentIndex(1)
                assert window.cards['Voyager'].values['speed'].text()=='Sin datos'
                assert not window.panel.markers['Voyager'].visible
                assert set(window.cards)=={'Hohmann','K','Voyager'}
                assert set(window.panel.positions)=={'Hohmann','K','Voyager','T','J','S'}
                window.seek(window.duration*.45);window.toggle();window.last_tick-=.1;window.tick();window.toggle()
                assert all(c.values['speed'].text() and c.values['saturn_distance'].text() for c in window.cards.values())
                for index in range(1,window.encounters.count()):
                    origin,target=window.encounters.itemData(index);window.encounters.setCurrentIndex(index)
                    assert abs(window.t-target)<1e-3 and window.panel.origin==origin
                    if origin=='T':
                        assert not window.panel.status.isHidden() and 'sin datos' in window.panel.status.text()
                        assert not np.any(window.panel.coords('T'))
                window.seek(window.duration*.45)
                assert all(button.isChecked() for button in window.visibility.values())
                for key,button in window.visibility.items():
                    button.click()
                    for visual in (window.panel.full,window.panel.trail,window.panel.markers,window.panel.labels):
                        assert not visual[key].visible
                        assert all(visual[other].visible for other in TRAJECTORIES if other!=key)
                    button.click()
                    assert window.panel.markers[key].visible
                for index in (1,2,3,0):window.frame.setCurrentIndex(index)
                window.panel.reset_camera();qt.processEvents()
                pixels=window.panel.canvas.render();assert np.std(pixels[:,:,:3])>3,'Render vacío'
                window.seek(window.duration);window.toggle();assert window.t==0 and window.playing;window.toggle()
                window.loop.setChecked(True);window.seek(window.duration-1);window.toggle();window.last_tick-=.1;window.tick()
                assert window.playing and window.t<window.duration;window.toggle()
                print('Escena superpuesta, tres trayectorias, tarjetas y reloj único comprobados.',flush=True)
            except Exception:
                import traceback
                traceback.print_exc();window.close();qt.exit(1);return
            window.close();qt.quit()
        QtCore.QTimer.singleShot(1500,verify)
    return qt.exec()
