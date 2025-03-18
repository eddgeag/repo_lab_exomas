#!/bin/bash

# Validación de parámetros
if [[ "$#" -lt 4 ]]; then
    echo "Uso: $0 -ref <ruta_referencia> -hilos <numero_hilos> -input_dir <ruta_de_entrada>"
    exit 1
fi

# Inicializar variables
REFERENCE=""
THREADS=""
MAIN_DIR=""

# Parsear parámetros
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -ref)
            REFERENCE="$2"
            shift 2
            ;;
        -hilos)
            THREADS="$2"
            shift 2
            ;;
        -input_dir)
            MAIN_DIR="$2"
            shift 2
            ;;
        *)
            echo "Parámetro desconocido: $1"
            exit 1
            ;;
    esac
done

# Validar parámetros obligatorios
if [[ -z "$REFERENCE"  || -z "$THREADS" || -z "$MAIN_DIR" ]]; then
    echo "Error: Faltan parámetros obligatorios."
    echo "Uso: $0 -ref <ruta_referencia>  -hilos <numero_hilos> -input_dir <ruta_de_entrada>"
    exit 1
fi

REFERENCE=$(realpath "$REFERENCE")
MAIN_DIR=$(realpath "$MAIN_DIR")
REFERENCE_DIR=$(dirname "$REFERENCE")

# Definir GATK_PATH para usar Docker
GATK_PATH_BASE="docker run --rm -v $(pwd):/data -v $REFERENCE_DIR:/reference -w /data"

# Validación de la referencia
if [[ ! -f "$REFERENCE" ]]; then
    echo "Archivo de referencia $REFERENCE no encontrado."
    exit 1
fi

# Crear diccionario de referencia si no existe
if [[ ! -f "${REFERENCE%.fa}.dict" ]]; then
    echo "Creando diccionario del fasta..."
    $GATK_PATH_BASE broadinstitute/gatk gatk CreateSequenceDictionary -R "/reference/$(basename $REFERENCE)"
fi

# Índice de referencia
if [[ -f "$REFERENCE.bwt" && -f "$REFERENCE.amb" && -f "$REFERENCE.ann" && -f "$REFERENCE.pac" && -f "$REFERENCE.sa" ]]; then
    echo "Índices de bwa ya existen. No se recrearán."
else
    echo "Creando índice bwa para el archivo de referencia..."
    bwa index "$REFERENCE"
fi

echo "Referencia: $REFERENCE"
echo "Índices: $REFERENCE.bwt, $REFERENCE.amb, $REFERENCE.ann, $REFERENCE.pac, $REFERENCE.sa"


perform_fastqc() {
    local input_dir="$1"   # Directorio donde están los archivos de entrada

    # Validar si se proporcionaron los parámetros obligatorios
    if [[ -z "$input_dir" ]]; then
        echo "Uso: perform_fastqc <input_dir>"
        return 1
    fi
     
    OUTPUT_DIR="$input_dir/pipeline_output"
    mkdir -p $OUTPUT_DIR
    
    # Buscar archivos con extensiones válidas
    local extensions=("*.fq.gz" "*.fastq.gz" "*.fq" "*.fastq")
    local files=()

    for ext in "${extensions[@]}"; do
        files+=($(find "$input_dir" -type f -name "$ext"))
    done

    # Verificar si se encontraron archivos
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No se encontraron archivos con extensiones válidas en $input_dir."
        return 1
    fi
    $output_dir_FQ="$OUTPUT_DIR/QC"
    mkdir $output_dir
    # Ejecutar FastQC
    echo "Realizando control de calidad con FastQC..."
    for file in "${files[@]}"; do
        echo "Procesando archivo: $file"
        fastqc -t 4 "$file" -o "$output_dir"
        if [[ $? -ne 0 ]]; then
            echo "Error al procesar $file con FastQC."
            return 1
        fi
    done

    echo "Control de calidad completado. Resultados guardados en $output_dir."
    return 0
}


