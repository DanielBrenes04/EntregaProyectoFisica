# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math


class SimuladorPlanoInclinado:
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador de segunda ley de newton con Friccion Estática")
        self.root.geometry("1200x800")

        # variables fisicas
        self.angulo = 30.0  # grados
        self.masa = 1.0  # kg
        self.coef_friccion = 0.15  # friccion cinetica base
        self.coef_friccion_estatica_base = 0.25  # friccion estática base
        self.fuerza_aplicada_magnitud = 0.0  # N
        self.fuerza_aplicada_angulo = 90.0  # grados respecto al plano
        self.g = 9.81  # m/s²

        # variables de la simulacion
        self.posicion_inicial_x = 1.0  # m - posicion inicial X
        self.posicion_inicial_y = 2.0  # m - posicion inicial Y (altura)
        self.posicion_x = 1.0  # m - posicion actual X
        self.posicion_y = 2.0  # m - posicion actual Y
        self.velocidad_x = 0.0  # m/s - velocidad en X
        self.velocidad_y = 0.0  # m/s - velocidad en Y
        self.aceleracion_x = 0.0  # m/s² - aceleracion en X
        self.aceleracion_y = 0.0  # m/s² - aceleracion en Y
        self.tiempo = 0.0  # s
        self.dt = 0.01
        self.simulando = False

        # Control de FPS
        self.fps_actual = 25  # FPS inicial (25 FPS = 40ms delay)
        self.fps_min = 5  # FPS mínimo
        self.fps_max = 60  # FPS máximo

        # Estados del bloque al inciio
        self.estado = "aire"

        # friccion
        self.friccion_estatica_activa = False
        self.coef_friccion_estatica_actual = self.coef_friccion_estatica_base

        # Configuracion del plano
        self.longitud_plano = 8.0  # m

        # Datos para graficos
        self.tiempo_datos = []
        self.posicion_x_datos = []
        self.posicion_y_datos = []
        self.velocidad_x_datos = []
        self.velocidad_y_datos = []
        self.velocidad_total_datos = []
        self.aceleracion_x_datos = []
        self.aceleracion_y_datos = []
        self.aceleracion_total_datos = []
        self.energia_datos = []
        self.friccion_datos = []  # Para graficar la friccion

        # Tiempo de inactividad iniciado en 0
        self.tiempo_sin_movimiento = 0.0

        self.crear_interfaz()
        self.determinar_estado()
        self.calcular_fuerzas()

    # UI
    def crear_interfaz(self):
        # Frame principal
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Panel de controles
        control_frame = ttk.LabelFrame(main_frame, text="Parametros de la Simulacion")
        control_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

        # Configuracion del angulo
        ttk.Label(control_frame, text="angulo de la pendiente (°):").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.angulo_var = tk.DoubleVar(value=self.angulo)
        angulo_scale = ttk.Scale(control_frame, from_=0, to=90, variable=self.angulo_var,
                                 command=self.actualizar_angulo, length=200)
        angulo_scale.grid(row=0, column=1, padx=5, pady=5)
        self.angulo_label = ttk.Label(control_frame, text=f"{self.angulo:.1f}°")
        self.angulo_label.grid(row=0, column=2, padx=5, pady=5)

        # Configuracion de la masa
        ttk.Label(control_frame, text="Masa del bloque (kg):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.masa_var = tk.DoubleVar(value=self.masa)
        masa_entry = ttk.Entry(control_frame, textvariable=self.masa_var, width=10)
        masa_entry.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        masa_entry.bind('<KeyRelease>', self.actualizar_masa)

        ttk.Label(control_frame, text="Coeficiente de friccion cinetica:").grid(row=2, column=0, sticky="w", padx=5,
                                                                                pady=5)
        self.friccion_var = tk.DoubleVar(value=self.coef_friccion)
        friccion_scale = ttk.Scale(control_frame, from_=0, to=0.30, variable=self.friccion_var,
                                   command=self.actualizar_friccion, length=200)
        friccion_scale.grid(row=2, column=1, padx=5, pady=5)
        self.friccion_label = ttk.Label(control_frame, text=f"{self.coef_friccion:.2f}")
        self.friccion_label.grid(row=2, column=2, padx=5, pady=5)

        ttk.Label(control_frame, text="Coeficiente de friccion estática base:").grid(row=3, column=0, sticky="w",
                                                                                     padx=5, pady=5)
        self.friccion_estatica_var = tk.DoubleVar(value=self.coef_friccion_estatica_base)
        friccion_estatica_scale = ttk.Scale(control_frame, from_=0, to=0.30, variable=self.friccion_estatica_var,
                                            command=self.actualizar_friccion_estatica, length=200)
        friccion_estatica_scale.grid(row=3, column=1, padx=5, pady=5)
        self.friccion_estatica_label = ttk.Label(control_frame, text=f"{self.coef_friccion_estatica_base:.2f}")
        self.friccion_estatica_label.grid(row=3, column=2, padx=5, pady=5)

        # Configuracion de fuerza aplicada
        ttk.Label(control_frame, text="Fuerza aplicada (N):").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.fuerza_mag_var = tk.DoubleVar(value=self.fuerza_aplicada_magnitud)
        fuerza_scale = ttk.Scale(control_frame, from_=0, to=100, variable=self.fuerza_mag_var,
                                 command=self.actualizar_fuerza, length=200)
        fuerza_scale.grid(row=4, column=1, padx=5, pady=5)
        self.fuerza_mag_label = ttk.Label(control_frame, text=f"{self.fuerza_aplicada_magnitud:.1f}")
        self.fuerza_mag_label.grid(row=4, column=2, padx=5, pady=5)

        ttk.Label(control_frame, text="angulo de fuerza aplicada (°):").grid(row=5, column=0, sticky="w", padx=5,
                                                                             pady=5)
        self.fuerza_ang_var = tk.DoubleVar(value=self.fuerza_aplicada_angulo)
        fuerza_ang_scale = ttk.Scale(control_frame, from_=0, to=180, variable=self.fuerza_ang_var,
                                     command=self.actualizar_angulo_fuerza, length=200)
        fuerza_ang_scale.grid(row=5, column=1, padx=5, pady=5)
        self.fuerza_ang_label = ttk.Label(control_frame, text=f"{self.fuerza_aplicada_angulo:.1f}°")
        self.fuerza_ang_label.grid(row=5, column=2, padx=5, pady=5)

        # Configuracion de posicion inicial
        ttk.Label(control_frame, text="Posicion inicial X (m):").grid(row=6, column=0, sticky="w", padx=5, pady=5)
        self.pos_inicial_x_var = tk.DoubleVar(value=self.posicion_inicial_x)
        pos_inicial_x_spinbox = tk.Spinbox(control_frame, from_=-1000.0, to=1000.0, increment=1.0,
                                           textvariable=self.pos_inicial_x_var, width=10,
                                           command=self.actualizar_posicion_inicial)
        pos_inicial_x_spinbox.grid(row=6, column=1, sticky="w", padx=5, pady=5)
        pos_inicial_x_spinbox.bind('<FocusOut>', self.actualizar_posicion_inicial)

        ttk.Label(control_frame, text="Posicion inicial Y (m):").grid(row=7, column=0, sticky="w", padx=5, pady=5)
        self.pos_inicial_y_var = tk.DoubleVar(value=self.posicion_inicial_y)
        pos_inicial_y_spinbox = tk.Spinbox(control_frame, from_=0.0, to=1000.0, increment=1.0,
                                           textvariable=self.pos_inicial_y_var, width=10,
                                           command=self.actualizar_posicion_inicial)
        pos_inicial_y_spinbox.grid(row=7, column=1, sticky="w", padx=5, pady=5)
        pos_inicial_y_spinbox.bind('<FocusOut>', self.actualizar_posicion_inicial)

        # Definir de longitud del plano
        ttk.Label(control_frame, text="Longitud del plano (m):").grid(row=8, column=0, sticky="w", padx=5, pady=5)
        self.longitud_var = tk.DoubleVar(value=self.longitud_plano)
        longitud_spinbox = tk.Spinbox(control_frame, from_=0.5, to=1000.0, increment=0.1,
                                      textvariable=self.longitud_var, width=10,
                                      command=self.actualizar_longitud_plano)
        longitud_spinbox.grid(row=8, column=1, sticky="w", padx=5, pady=5)
        longitud_spinbox.bind('<FocusOut>', self.actualizar_longitud_plano)

        # Botones
        button_frame = ttk.Frame(control_frame)
        button_frame.grid(row=9, column=0, columnspan=3, pady=10)

        self.start_button = ttk.Button(button_frame, text="Iniciar/Resumir simulacion", command=self.iniciar_simulacion)
        self.start_button.pack(side=tk.LEFT, padx=5)

        self.stop_button = ttk.Button(button_frame, text="Detener", command=self.detener_simulacion)
        self.stop_button.pack(side=tk.LEFT, padx=5)

        self.reset_button = ttk.Button(button_frame, text="Reiniciar", command=self.reiniciar_simulacion)
        self.reset_button.pack(side=tk.LEFT, padx=5)

        self.default_button = ttk.Button(button_frame, text="Valores predeterminados",
                                         command=self.restablecer_valores_predeterminados)
        self.default_button.pack(side=tk.LEFT, padx=5)

        speed_frame = ttk.Frame(control_frame)
        speed_frame.grid(row=10, column=0, columnspan=3, pady=5)

        ttk.Label(speed_frame, text="Velocidad de simulacion:").pack(side=tk.LEFT, padx=5)

        self.fps_down_button = ttk.Button(speed_frame, text="←", command=self.disminuir_fps, width=8)
        self.fps_down_button.pack(side=tk.LEFT, padx=2)

        self.fps_label = ttk.Label(speed_frame, text=f"{self.fps_actual} FPS", width=10)
        self.fps_label.pack(side=tk.LEFT, padx=5)

        self.fps_up_button = ttk.Button(speed_frame, text="→", command=self.aumentar_fps, width=8)
        self.fps_up_button.pack(side=tk.LEFT, padx=2)

        # Panel de informacion
        info_frame = ttk.LabelFrame(control_frame, text="Estado Actual")
        info_frame.grid(row=11, column=0, columnspan=3, sticky="ew", pady=10)

        self.info_text = tk.Text(info_frame, height=35, width=40)
        self.info_text.pack(padx=5, pady=5)

        # Frame para graficos
        graph_frame = ttk.Frame(main_frame)
        graph_frame.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)

        # Configurar el grafico
        self.fig = plt.figure(figsize=(8, 8))
        self.ax1 = self.fig.add_subplot(2, 1, 1, aspect='equal')
        self.ax2 = self.fig.add_subplot(2, 2, 3)
        self.ax3 = self.fig.add_subplot(2, 2, 4)

        self.canvas = FigureCanvasTkAgg(self.fig, graph_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)

        self.actualizar_graficos()

    # Logica
    def actualizar_angulo(self, value=None):
        self.angulo = float(value) if value else self.angulo_var.get()
        self.angulo_label.config(text=f"{self.angulo:.1f}°")
        self.determinar_estado()
        self.calcular_fuerzas()

    # Logica
    def actualizar_masa(self, event=None):
        try:
            self.masa = self.masa_var.get()
            self.calcular_fuerzas()
        except:
            pass

    # Logica
    def actualizar_friccion(self, value=None):
        self.coef_friccion = float(value) if value else self.friccion_var.get()
        self.friccion_label.config(text=f"{self.coef_friccion:.2f}")
        self.calcular_fuerzas()

    # Logica para friccion estática
    def actualizar_friccion_estatica(self, value=None):
        self.coef_friccion_estatica_base = float(value) if value else self.friccion_estatica_var.get()
        self.friccion_estatica_label.config(text=f"{self.coef_friccion_estatica_base:.2f}")
        self.calcular_fuerzas()

    # Logica
    def actualizar_fuerza(self, value=None):
        try:
            self.fuerza_aplicada_magnitud = float(value) if value else self.fuerza_mag_var.get()
            self.fuerza_mag_label.config(text=f"{self.fuerza_aplicada_magnitud:.1f}")
            self.calcular_fuerzas()
        except:
            pass

    # Logica
    def actualizar_angulo_fuerza(self, value=None):
        self.fuerza_aplicada_angulo = float(value) if value else self.fuerza_ang_var.get()
        self.fuerza_ang_label.config(text=f"{self.fuerza_aplicada_angulo:.1f}°")
        self.calcular_fuerzas()

    # Logica
    def actualizar_posicion_inicial(self, event=None):
        try:
            nueva_pos_x = self.pos_inicial_x_var.get()
            nueva_pos_y = self.pos_inicial_y_var.get()

            if nueva_pos_y < 0:
                messagebox.showwarning("Advertencia", "La altura inicial debe ser mayor o igual a 0")
                self.pos_inicial_y_var.set(self.posicion_inicial_y)
                return

            self.posicion_inicial_x = nueva_pos_x
            self.posicion_inicial_y = nueva_pos_y

            # Guardar posicion
            if not self.simulando:
                self.posicion_x = self.posicion_inicial_x
                self.posicion_y = self.posicion_inicial_y
                self.determinar_estado()
            self.calcular_fuerzas()
        except:
            pass

    # Logica
    def actualizar_longitud_plano(self, event=None):
        try:
            nueva_longitud = self.longitud_var.get()
            if nueva_longitud > 0:
                self.longitud_plano = nueva_longitud
                self.determinar_estado()
                self.calcular_fuerzas()
            else:
                messagebox.showwarning("Advertencia", "La longitud del plano debe ser mayor a 0")
                self.longitud_var.set(self.longitud_plano)
        except:
            pass

    # Logica
    def determinar_estado(self):
        estado_anterior = self.estado

        if self.posicion_x <= 0:
            altura_plano_en_x = 0
        elif self.posicion_x <= self.longitud_plano:
            altura_plano_en_x = self.posicion_x * math.tan(math.radians(self.angulo))
        else:
            altura_plano_en_x = 0  # y=0 es suelo

        tolerancia = 0.05  # Tolerancia de cambio de friccion

        if self.posicion_y > altura_plano_en_x + tolerancia:
            self.estado = "aire"
        elif (self.posicion_x >= -tolerancia and self.posicion_x <= self.longitud_plano + tolerancia and
              abs(self.posicion_y - altura_plano_en_x) <= tolerancia):
            self.estado = "plano"
        elif abs(self.posicion_y) <= tolerancia:
            self.estado = "suelo"
        else:
            self.estado = "aire"

        # TEsting: Imprimir cuando hay cambio de estado
        if estado_anterior != self.estado:
            print(
                f"Cambio de estado: {estado_anterior} -> {self.estado} en t={self.tiempo:.3f}s, pos=({self.posicion_x:.3f}, {self.posicion_y:.3f})")
            print(f"  Velocidad: ({self.velocidad_x:.3f}, {self.velocidad_y:.3f})")
            print(f"  Altura plano en X: {altura_plano_en_x:.3f}")
            print("---")

    def es_subiendo_pendiente(self):
        if self.estado != "plano":
            return False

        angulo_rad = math.radians(self.angulo)
        # Velocidad  plano
        velocidad_paralela = self.velocidad_x * math.cos(angulo_rad) + self.velocidad_y * math.sin(angulo_rad)

        return velocidad_paralela > 0.01  # Umbral pequeño para evitar ruido numerico

    def calcular_coeficiente_friccion_estatica_dinamico(self, fuerza_aplicada_magnitud):
        # Incrementar segun aumenta la fuerza
        factor_incremento = 0.05
        incremento = fuerza_aplicada_magnitud * factor_incremento
        coef_max = self.coef_friccion_estatica_base * 2.0

        coef_dinamico = min(self.coef_friccion_estatica_base + incremento, coef_max)
        return coef_dinamico

    # Logica
    def calcular_fuerzas(self):
        # Calcular fuerzas y aceleracion del bloque
        self.aceleracion_x = 0.0
        self.aceleracion_y = 0.0
        friccion_aplicada = 0.0

        if self.estado == "aire":
            # Caida libre
            self.aceleracion_x = 0.0
            self.aceleracion_y = -self.g
            self.friccion_estatica_activa = False

        elif self.estado == "plano":
            # Plano inclinado
            angulo_rad = math.radians(self.angulo)
            fuerza_ang_rad = math.radians(self.fuerza_aplicada_angulo)

            # peso
            peso_paralelo = self.masa * self.g * math.sin(angulo_rad)
            peso_perpendicular = self.masa * self.g * math.cos(angulo_rad)

            # fuerza aplicada en plano
            fuerza_paralelo = self.fuerza_aplicada_magnitud * math.cos(fuerza_ang_rad)
            fuerza_perpendicular = self.fuerza_aplicada_magnitud * math.sin(fuerza_ang_rad)

            # Fuerza normal
            normal = peso_perpendicular + fuerza_perpendicular

            # Determinar si está subiendo la pendiente
            subiendo = self.es_subiendo_pendiente()

            if subiendo:
                # Friccion estatica
                self.coef_friccion_estatica_actual = self.calcular_coeficiente_friccion_estatica_dinamico(
                    self.fuerza_aplicada_magnitud)
                friccion_maxima_estatica = self.coef_friccion_estatica_actual * abs(normal)

                fuerza_neta_sin_friccion = fuerza_paralelo - peso_paralelo

                # Aplicar friccion estática
                if abs(fuerza_neta_sin_friccion) <= friccion_maxima_estatica and abs(self.velocidad_x) < 0.01:
                    aceleracion_paralela = 0.0
                    friccion_aplicada = -fuerza_neta_sin_friccion  # Friccion exacta para equilibrar
                    self.friccion_estatica_activa = True
                    print(
                        f"Friccion estática pura: {friccion_aplicada:.3f} N, μₛ = {self.coef_friccion_estatica_actual:.3f}")
                else:
                    # Friccion estática dinámica
                    friccion_real = -friccion_maxima_estatica
                    fuerza_neta = fuerza_neta_sin_friccion + friccion_real
                    aceleracion_paralela = fuerza_neta / self.masa
                    friccion_aplicada = friccion_real
                    self.friccion_estatica_activa = True
                    print(
                        f"Friccion estática dinámica: {friccion_aplicada:.3f} N, μₛ = {self.coef_friccion_estatica_actual:.3f}")

            else:
                # Friccion cinetica
                friccion_maxima = self.coef_friccion * abs(normal)
                fuerza_neta_sin_friccion = fuerza_paralelo - peso_paralelo

                if abs(fuerza_neta_sin_friccion) <= friccion_maxima and abs(self.velocidad_x) < 0.01:
                    aceleracion_paralela = 0.0
                    friccion_aplicada = -fuerza_neta_sin_friccion
                else:
                    if fuerza_neta_sin_friccion > 0:
                        friccion_real = -friccion_maxima
                    else:
                        friccion_real = friccion_maxima

                    fuerza_neta = fuerza_neta_sin_friccion + friccion_real
                    aceleracion_paralela = fuerza_neta / self.masa
                    friccion_aplicada = friccion_real

                self.friccion_estatica_activa = False

            self.aceleracion_x = aceleracion_paralela * math.cos(angulo_rad)
            self.aceleracion_y = aceleracion_paralela * math.sin(angulo_rad)

        elif self.estado == "suelo":
            fuerza_ang_rad = math.radians(self.fuerza_aplicada_angulo)

            peso_perpendicular = self.masa * self.g
            fuerza_horizontal = self.fuerza_aplicada_magnitud * math.cos(fuerza_ang_rad)
            fuerza_vertical = self.fuerza_aplicada_magnitud * math.sin(fuerza_ang_rad)

            normal = peso_perpendicular + fuerza_vertical
            friccion_maxima = self.coef_friccion * abs(normal)

            if abs(self.velocidad_x) > 0.01:
                friccion_real = -friccion_maxima if self.velocidad_x > 0 else friccion_maxima
                fuerza_neta_x = fuerza_horizontal + friccion_real
                self.aceleracion_x = fuerza_neta_x / self.masa
                self.aceleracion_x = self.aceleracion_x * 0.6
                friccion_aplicada = friccion_real
            else:  # Sin movimiento significativo
                if abs(fuerza_horizontal) <= friccion_maxima:
                    self.aceleracion_x = 0.0
                    friccion_aplicada = -fuerza_horizontal
                else:
                    friccion_real = -friccion_maxima if fuerza_horizontal > 0 else friccion_maxima
                    fuerza_neta_x = fuerza_horizontal + friccion_real
                    self.aceleracion_x = fuerza_neta_x / self.masa
                    friccion_aplicada = friccion_real

            self.aceleracion_y = 0.0
            self.friccion_estatica_activa = False

        # datos de friccion para gráficos
        if hasattr(self, 'friccion_datos'):
            self.friccion_datos.append(friccion_aplicada)

        self.actualizar_info()

    # UI
    def actualizar_info(self):
        energia_cinetica = 0.5 * self.masa * (self.velocidad_x ** 2 + self.velocidad_y ** 2)
        energia_potencial = self.masa * self.g * self.posicion_y
        energia_total = energia_cinetica + energia_potencial
        velocidad_total = math.sqrt(self.velocidad_x ** 2 + self.velocidad_y ** 2)

        # Informacion adicional sobre friccion
        subiendo = self.es_subiendo_pendiente() if self.estado == "plano" else False
        tipo_friccion = "Estática" if self.friccion_estatica_activa else "Cinetica"

        info = f"""ESTADO DEL BLOQUE:
        Estado actual: {self.estado.upper()}
        {"↑ SUBIENDO PENDIENTE" if subiendo else "↓ Bajando/Reposo" if self.estado == "plano" else ""}

        POSICION:
        X: {self.posicion_x:.2f} m
        Y: {self.posicion_y:.2f} m

        VELOCIDAD:
        Vx: {abs(self.velocidad_x):.2f} m/s
        Vy: {abs(self.velocidad_y):.2f} m/s
        V total: {velocidad_total:.2f} m/s

        ACELERACION:
        Ax: {abs(self.aceleracion_x):.2f} m/s²
        Ay: {abs(self.aceleracion_y):.2f} m/s²

        FRICCION:
        Tipo: {tipo_friccion}
        μ cinetico: {self.coef_friccion:.3f}
        μ estático base: {self.coef_friccion_estatica_base:.3f}
        {"μ estático actual: " + f"{self.coef_friccion_estatica_actual:.3f}" if self.friccion_estatica_activa else ""}

        ENERGIA:
        Cinetica: {energia_cinetica:.2f} J
        Potencial: {energia_potencial:.2f} J
        Total: {energia_total:.2f} J

        CONFIGURACION:
        Masa: {self.masa:.1f} kg
        Tiempo: {self.tiempo:.2f} s
                """

        self.info_text.delete(1.0, tk.END)
        self.info_text.insert(1.0, info)

    # UI
    def iniciar_simulacion(self):
        # Start
        if not self.simulando:
            self.simulando = True
            self.start_button.config(state="disabled")
            self.animar()

    # UI
    def detener_simulacion(self):
        # Stop
        self.simulando = False
        self.start_button.config(state="normal")

    # UI
    def reiniciar_simulacion(self):
        # Reset
        self.detener_simulacion()
        self.posicion_x = self.posicion_inicial_x
        self.posicion_y = self.posicion_inicial_y
        self.velocidad_x = 0.0
        self.velocidad_y = 0.0
        self.tiempo = 0.0
        self.tiempo_datos = []
        self.posicion_x_datos = []
        self.posicion_y_datos = []
        self.velocidad_x_datos = []
        self.velocidad_y_datos = []
        self.velocidad_total_datos = []
        self.aceleracion_x_datos = []
        self.aceleracion_y_datos = []
        self.aceleracion_total_datos = []
        self.energia_datos = []
        self.friccion_datos = []
        self.tiempo_sin_movimiento = 0.0
        self.friccion_estatica_activa = False
        self.determinar_estado()
        self.calcular_fuerzas()
        self.actualizar_graficos()

    def restablecer_valores_predeterminados(self):
        # Valores iniciales predeterminados
        self.angulo = 30.0
        self.masa = 1.0
        self.coef_friccion = 0.15
        self.coef_friccion_estatica_base = 0.25
        self.fuerza_aplicada_magnitud = 0.0
        self.fuerza_aplicada_angulo = 90.0
        self.posicion_inicial_x = 1.0
        self.posicion_inicial_y = 2.0
        self.longitud_plano = 8.0

        self.angulo_var.set(self.angulo)
        self.angulo_label.config(text=f"{self.angulo:.1f}°")

        self.masa_var.set(self.masa)

        self.friccion_var.set(self.coef_friccion)
        self.friccion_label.config(text=f"{self.coef_friccion:.2f}")

        self.friccion_estatica_var.set(self.coef_friccion_estatica_base)
        self.friccion_estatica_label.config(text=f"{self.coef_friccion_estatica_base:.2f}")

        self.fuerza_mag_var.set(self.fuerza_aplicada_magnitud)

        self.fuerza_ang_var.set(self.fuerza_aplicada_angulo)
        self.fuerza_ang_label.config(text=f"{self.fuerza_aplicada_angulo:.1f}°")

        self.pos_inicial_x_var.set(self.posicion_inicial_x)
        self.pos_inicial_y_var.set(self.posicion_inicial_y)

        self.longitud_var.set(self.longitud_plano)

        self.fps_actual = 25
        self.fps_label.config(text=f"{self.fps_actual} FPS")

        # Reiniciar a valores predeterminados
        self.reiniciar_simulacion()

    def disminuir_fps(self):
        if self.fps_actual > self.fps_min:
            self.fps_actual = max(self.fps_min, self.fps_actual - 5)
            self.fps_label.config(text=f"{self.fps_actual} FPS")
            print(f"FPS reducidos a: {self.fps_actual}")

    def aumentar_fps(self):
        if self.fps_actual < self.fps_max:
            self.fps_actual = min(self.fps_max, self.fps_actual + 5)
            self.fps_label.config(text=f"{self.fps_actual} FPS")
            print(f"FPS aumentados a: {self.fps_actual}")

    def calcular_delay_animacion(self):
        return int(1000 / self.fps_actual)

    # UI
    def paso_simulacion(self):
        if self.simulando:
            self.calcular_fuerzas()

            nueva_velocidad_x = self.velocidad_x + self.aceleracion_x * self.dt
            nueva_velocidad_y = self.velocidad_y + self.aceleracion_y * self.dt

            nueva_posicion_x = self.posicion_x + nueva_velocidad_x * self.dt
            nueva_posicion_y = self.posicion_y + nueva_velocidad_y * self.dt

            self.verificar_colisiones(nueva_posicion_x, nueva_posicion_y, nueva_velocidad_x, nueva_velocidad_y)

            # coreccion para que no atraviese el suelo
            if self.estado == "suelo" and self.posicion_y < 0:
                self.posicion_y = 0
                self.velocidad_y = 0

            self.tiempo += self.dt

            # Guardar datos para graficos
            self.tiempo_datos.append(self.tiempo)
            self.posicion_x_datos.append(self.posicion_x)
            self.posicion_y_datos.append(self.posicion_y)
            self.velocidad_x_datos.append(self.velocidad_x)
            self.velocidad_y_datos.append(self.velocidad_y)

            velocidad_total = math.sqrt(self.velocidad_x ** 2 + self.velocidad_y ** 2)
            self.velocidad_total_datos.append(velocidad_total)

            # Guardar datos de aceleracion
            self.aceleracion_x_datos.append(self.aceleracion_x)
            self.aceleracion_y_datos.append(self.aceleracion_y)
            aceleracion_total = math.sqrt(self.aceleracion_x ** 2 + self.aceleracion_y ** 2)
            self.aceleracion_total_datos.append(aceleracion_total)

            energia_total = 0.5 * self.masa * (
                    self.velocidad_x ** 2 + self.velocidad_y ** 2) + self.masa * self.g * self.posicion_y
            self.energia_datos.append(energia_total)

            # Detectar inactividad
            if velocidad_total < 1:
                self.tiempo_sin_movimiento += self.dt
            else:
                self.tiempo_sin_movimiento = 0.0

            if self.tiempo_sin_movimiento >= 1.0:  # 1 segundo sin movimiento
                self.detener_simulacion()
                print(f"Simulacion detenida por inactividad. Velocidad: {velocidad_total:.3f} m/s")

            # Limitar datos para eviat consumo de memoria/crash
            if len(self.tiempo_datos) > 2000:
                self.tiempo_datos.pop(0)
                self.posicion_x_datos.pop(0)
                self.posicion_y_datos.pop(0)
                self.velocidad_x_datos.pop(0)
                self.velocidad_y_datos.pop(0)
                self.velocidad_total_datos.pop(0)
                self.aceleracion_x_datos.pop(0)
                self.aceleracion_y_datos.pop(0)
                self.aceleracion_total_datos.pop(0)
                self.energia_datos.pop(0)
                if len(self.friccion_datos) > 0:
                    self.friccion_datos.pop(0)

    # Logica
    def verificar_colisiones(self, nueva_x, nueva_y, nueva_vx, nueva_vy):
        # Calcular altura del plano en la nueva posicion X
        if nueva_x <= 0:
            altura_plano = 0
        elif nueva_x <= self.longitud_plano:
            altura_plano = nueva_x * math.tan(math.radians(self.angulo))
        else:
            altura_plano = 0

        # Verificar transicion al suelo
        if nueva_y <= 0:
            self.posicion_x = nueva_x
            self.posicion_y = 0  # Mantenerse en suelo
            self.velocidad_x = nueva_vx  # Conservar velocidad horizontal
            self.velocidad_y = 0  # Anular velocidad vertical
            self.estado = "suelo"
            return

        # validar colision con el plano
        if (nueva_x >= 0 and nueva_x <= self.longitud_plano and
                nueva_y <= altura_plano and self.posicion_y > altura_plano):
            angulo_rad = math.radians(self.angulo)

            # Componente de velocidad paralela en plano
            vel_paralela = (nueva_vx * math.cos(angulo_rad) + nueva_vy * math.sin(angulo_rad))

            # Componente perpendicular se anula en la colision
            self.velocidad_x = vel_paralela * math.cos(angulo_rad)
            self.velocidad_y = vel_paralela * math.sin(angulo_rad)

            self.posicion_x = nueva_x
            self.posicion_y = altura_plano
            self.estado = "plano"
            return

        # actualizar normalmente
        self.posicion_x = nueva_x
        self.posicion_y = nueva_y
        self.velocidad_x = nueva_vx
        self.velocidad_y = nueva_vy

        self.determinar_estado()

    # calcular componentes de fuerza aplicada
    def calcular_componentes_fuerza_global(self):
        if self.fuerza_aplicada_magnitud == 0:
            return 0, 0

        if self.estado == "plano":
            # Cuando está en el plano, el ángulo de fuerza es respecto al plano inclinado
            angulo_plano_rad = math.radians(self.angulo)
            fuerza_ang_rad = math.radians(self.fuerza_aplicada_angulo)

            # angulo total respecto al eje X global
            angulo_total = angulo_plano_rad + fuerza_ang_rad

            fx_global = self.fuerza_aplicada_magnitud * math.cos(angulo_total)
            fy_global = self.fuerza_aplicada_magnitud * math.sin(angulo_total)
        else:
            # En aire o suelo, el angulo de fuerza es respecto al eje X horizontal
            fuerza_ang_rad = math.radians(self.fuerza_aplicada_angulo)
            fx_global = self.fuerza_aplicada_magnitud * math.cos(fuerza_ang_rad)
            fy_global = self.fuerza_aplicada_magnitud * math.sin(fuerza_ang_rad)

        return fx_global, fy_global

    # UI
    def actualizar_graficos(self):
        # Reiniciar graficos
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()

        # Grafico del plano inclinado
        margen_vista = 4
        # Rango fijo para la camara
        rango = 2 * margen_vista  # ventana de 8 unidades

        # Centrar la vista en la posicion actual del bloque
        x_centro = self.posicion_x
        y_centro = self.posicion_y

        # Limites cuadrados para evitar deformacion
        self.ax1.set_xlim(x_centro - margen_vista, x_centro + margen_vista)
        self.ax1.set_ylim(y_centro - margen_vista, y_centro + margen_vista)
        self.ax1.set_aspect('equal')

        # Dibujar lineas de los ejes X y Y
        self.ax1.axhline(0, color='black', linewidth=1)

        # Dibujar suelo
        max_x = max(self.longitud_plano + margen_vista, self.posicion_x + margen_vista)
        x_suelo = [-margen_vista, max_x]
        y_suelo = [0, 0]
        self.ax1.axhline(y=0, color='k', linewidth=3, label='Superficie horizontal (suelo)')

        # Dibujar plano inclinado
        x_plano = [0, self.longitud_plano]
        y_plano = [0, self.longitud_plano * math.tan(math.radians(self.angulo))]
        self.ax1.plot(x_plano, y_plano, 'k-', linewidth=3, label='Plano inclinado')

        # color de friccion estatica
        if self.friccion_estatica_activa and self.es_subiendo_pendiente():
            color_bloque = 'orange'  # Color especial para friccion estática
            label_bloque = f'Bloque ({self.estado}) - FRICCION ESTATICA'
        else:
            color_bloque = {'aire': 'red', 'plano': 'blue', 'suelo': 'green'}[self.estado]
            label_bloque = f'Bloque ({self.estado})'

        self.ax1.plot(self.posicion_x, self.posicion_y, 's', color=color_bloque,
                      markersize=15, label=label_bloque)

        # Dibujar vector de fuerza aplicada
        fx_global, fy_global = self.calcular_componentes_fuerza_global()

        if self.fuerza_aplicada_magnitud > 0.1:
            # Escalar la flecha para visualizacion
            escala_fuerza = 0.8
            fx_display = fx_global * escala_fuerza
            fy_display = fy_global * escala_fuerza

            # Dibujar la flecha de fuerza aplicada
            self.ax1.arrow(self.posicion_x, self.posicion_y, fx_display, fy_display,
                           head_width=0.25, head_length=0.25, fc='orange', ec='darkorange',
                           linewidth=4, alpha=0.8, label=f'Fuerza aplicada ({self.fuerza_aplicada_magnitud:.1f} N)')

        # Dibujar vectores de velocidad
        if len(self.tiempo_datos) > 0:
            velocidad_total = math.sqrt(self.velocidad_x ** 2 + self.velocidad_y ** 2)
            if velocidad_total > 0.1:
                escala = 0.5
                vx_display = self.velocidad_x * escala
                vy_display = self.velocidad_y * escala
                self.ax1.arrow(self.posicion_x, self.posicion_y, vx_display, vy_display,
                               head_width=0.2, head_length=0.2, fc='purple', ec='purple',
                               label=f'Velocidad ({velocidad_total:.1f} m/s)')

        # trayectoria
        if len(self.posicion_x_datos) > 1:
            self.ax1.plot(self.posicion_x_datos, self.posicion_y_datos, 'r--', alpha=0.5, label='Trayectoria')

        self.ax1.set_xlabel('Posicion X (m)')
        self.ax1.set_ylabel('Posicion Y (m)')
        self.ax1.set_title('Simulacion del Plano Inclinado con Friccion Estatica')
        self.ax1.legend()
        self.ax1.grid(True, alpha=0.3)

        if len(self.tiempo_datos) > 1:
            margen_tiempo = 0.5
            # Grafico 1: Magnitudes de velocidad y aceleracion en X
            velocidad_x_magnitud = [abs(v) for v in self.velocidad_x_datos]
            aceleracion_x_magnitud = [abs(a) for a in self.aceleracion_x_datos]

            self.ax2.plot(self.tiempo_datos, velocidad_x_magnitud, 'r-', label='Velocidad X (magnitud, m/s)')
            self.ax2.plot(self.tiempo_datos, aceleracion_x_magnitud, 'b-', label='Aceleracion X (magnitud, m/s²)')

            # Friccion estatica esta activa
            if len(self.friccion_datos) == len(self.tiempo_datos):
                friccion_estatica_momentos = []
                for i, t in enumerate(self.tiempo_datos):
                    # Determinar si habia friccion estatica en ese momento
                    if i < len(self.posicion_x_datos):
                        x_pos = self.posicion_x_datos[i]
                        # Si estaba en el plano y habia friccion significativa
                        if (x_pos > 0 and x_pos < self.longitud_plano and
                                len(self.friccion_datos) > i and abs(self.friccion_datos[i]) > 0.1):
                            friccion_estatica_momentos.append(t)

                # Sombreo de friccion estatica
                for t in friccion_estatica_momentos:
                    self.ax2.axvline(x=t, color='orange', alpha=0.3, linewidth=0.5)

            self.ax2.plot(self.tiempo_datos[-1], velocidad_x_magnitud[-1], 'ro', markersize=10)
            self.ax2.plot(self.tiempo_datos[-1], aceleracion_x_magnitud[-1], 'bo', markersize=10)
            self.ax2.set_xlabel('Tiempo (s)')
            self.ax2.set_ylabel('Magnitud')
            self.ax2.set_title('Velocidad y Aceleracion en X (Magnitudes)')
            self.ax2.legend()
            self.ax2.grid(True, alpha=0.3)
            self.ax2.set_xlim(min(self.tiempo_datos), max(self.tiempo_datos) + margen_tiempo)

            # Grafico 2: Magnitudes de velocidad y aceleracion en Y
            velocidad_y_magnitud = [abs(v) for v in self.velocidad_y_datos]
            aceleracion_y_magnitud = [abs(a) for a in self.aceleracion_y_datos]

            self.ax3.plot(self.tiempo_datos, velocidad_y_magnitud, 'r-', label='Velocidad Y (magnitud, m/s)')
            self.ax3.plot(self.tiempo_datos, aceleracion_y_magnitud, 'b-', label='Aceleracion Y (magnitud, m/s²)')
            self.ax3.plot(self.tiempo_datos[-1], velocidad_y_magnitud[-1], 'ro', markersize=10)
            self.ax3.plot(self.tiempo_datos[-1], aceleracion_y_magnitud[-1], 'bo', markersize=10)

            # Linea horizontal en cero para referencia
            self.ax3.axhline(y=0, color='black', linestyle='-', alpha=0.3)

            self.ax3.set_xlabel('Tiempo (s)')
            self.ax3.set_ylabel('Magnitud')
            self.ax3.set_title('Velocidad y Aceleracion en Y (Magnitudes)')
            self.ax3.legend()
            self.ax3.grid(True, alpha=0.3)
            self.ax3.set_xlim(min(self.tiempo_datos), max(self.tiempo_datos) + margen_tiempo)

        # Deteccion de estado
        estado_anterior = None
        for i, t in enumerate(self.tiempo_datos):
            if i < len(self.posicion_x_datos):
                x_pos = self.posicion_x_datos[i]
                y_pos = self.posicion_y_datos[i]

                if y_pos > 0.1:
                    estado_actual = "aire"
                elif x_pos <= self.longitud_plano and abs(y_pos - x_pos * math.tan(math.radians(self.angulo))) < 0.1:
                    estado_actual = "plano"
                else:
                    estado_actual = "suelo"

                if estado_anterior is not None and estado_actual != estado_anterior:
                    self.ax2.axvline(x=t, color='gray', linestyle='--', alpha=0.7)
                    self.ax3.axvline(x=t, color='gray', linestyle='--', alpha=0.7)

                estado_anterior = estado_actual

        self.canvas.draw()

    # UI
    def animar(self):
        if self.simulando:
            self.paso_simulacion()
            self.actualizar_graficos()

            # Usar FPS dinamicos
            delay = self.calcular_delay_animacion()
            self.root.after(delay, self.animar)


def main():
    root = tk.Tk()
    app = SimuladorPlanoInclinado(root)
    root.mainloop()


if __name__ == "__main__":
    main()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
