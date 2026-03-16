---
description: How to execute Google Colab notebooks from Antigravity IDE
---

# Google Colab Integration in Antigravity

This workflow guides you through setting up the official Google Colab extension to run notebooks directly from Antigravity.

## Step 1: Install Required Extensions

1. Open the **Extensions** view in Antigravity (`Ctrl+Shift+X`).
2. Search for and install the following extensions:
   - **Jupyter** (by Microsoft) - *Required dependency*
   - **Google Colab** (by Google)
3. Restart or reload Antigravity if prompted.

## Step 2: Create or Open a Notebook

1. Navigate to your project folder.
2. Create a new file with the `.ipynb` extension (e.g., `test_colab.ipynb`) or open an existing one.

## Step 3: Connect to Google Colab Runtime

1. Open your `.ipynb` file.
2. Click on **Select Kernel** in the top-right corner of the editor.
3. Select **Select Another Kernel...** from the dropdown.
4. Choose **Google Colab**.
5. A browser window will open asking you to sign in with your Google account.
6. Once signed in, you will be prompted to select a Colab server:
   - Choose **New Colab Server** to start a fresh instance.
   - Choose your desired accelerator (CPU, GPU, or TPU) based on your Colab subscription.

## Step 4: Execute Code

1. Type some Python code in a cell (e.g., `print("Hello from Colab!")`).
2. Run the cell (`Shift+Enter`).
3. The code will execute on the remote Google Colab infrastructure.

---

# Integración de Google Colab en Antigravity

Esta guía te ayuda a configurar la extensión oficial de Google Colab para ejecutar libretas directamente desde Antigravity.

## Paso 1: Instalar Extensiones Requeridas

1. Abre la vista de **Extensiones** en Antigravity (`Ctrl+Shift+X`).
2. Busca e instala las siguientes extensiones:
   - **Jupyter** (de Microsoft) - *Dependencia necesaria*
   - **Google Colab** (de Google)
3. Reinicia o recarga Antigravity si se te solicita.

## Paso 2: Crear o Abrir una Libreta

1. Navega a tu carpeta de proyecto.
2. Crea un nuevo archivo con extensión `.ipynb` (ej. `prueba_colab.ipynb`) o abre uno existente.

## Paso 3: Conectar al Runtime de Google Colab

1. Abre tu archivo `.ipynb`.
2. Haz clic en **Select Kernel** (Seleccionar Kernel) en la esquina superior derecha del editor.
3. Elige **Select Another Kernel...** (Seleccionar otro Kernel...) en el menú desplegable.
4. Selecciona **Google Colab**.
5. Se abrirá una ventana del navegador pidiéndote que inicies sesión con tu cuenta de Google.
6. Una vez que hayas iniciado sesión, se te pedirá que selecciones un servidor de Colab:
   - Elige **New Colab Server** (Nuevo servidor de Colab) para iniciar una instancia nueva.
   - Selecciona el acelerador deseado (CPU, GPU o TPU) según tu suscripción de Colab.

## Paso 4: Ejecutar Código

1. Escribe algo de código Python en una celda (ej. `print("¡Hola desde Colab!")`).
2. Ejecuta la celda (`Shift+Enter`).
3. El código se ejecutará en la infraestructura remota de Google Colab.
