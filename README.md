# csenergy-python
Library for Concentrated Solar Energy Systems


El usuario selecciona:
    
    - Fichero Weather, se crea dataframe con registros a intervalos temporales
    - Configuración de planta (N_campos, N_lazos, N_SCA_lazo, N_HCE_SCA)
    - Modelo de planta:
        - Modelo 1: Parámetros sin dispersión
        - Modelo 2: Parámetros con dispersión
    - Tipo de simulacion:
        - Simulacion 1: Comparación con datos de generacion
        - Simulacion 2: Simulación de planta
    
    - Se genera la planta:
        - Modelo 1: Por ser todos iguales no sería necesario crear N_lazosxN_SCA_lazoxN_HCE_SCA
        sino que se va a recurrir a bucles para generar las salidas de temperatura y guardarlos
        en el dataframe
        
        - Modelo 2: Se crean los N_lazosxN_SCA_lazoxN_HCE_SCA HCE con sus atributos cargados según 
        las funciones de dispersión que toque. 
        
    - Simulación:
        - Simulación 1: para cada registro temporal se toma el caudal real y se calcula la temperatura de
        salida con el modelo seleccionado. Hay que manejar las temperaturas por encima del setpoint como pérdidas de energía por 
        culpa del desenfoque
        
        - Simulación 2:  para cada registro temporal se calcula el caudal para conseguir la temperatura
        de salida.
        
        datos que debe almacenarse en cada HCE
        - Rendimiento, tin, tout, massflow
        
        datos que debe almacenarse en cada SCA
        - Rendimiento (promedio de sus HCE), tin, tout, t_probe,  massflow, status (desenfoque, etc...)
        
        datos que debe almacenarse en cada LOOP
        - Rendimiento (promedio de sus SCA), tin, tout, massflow
