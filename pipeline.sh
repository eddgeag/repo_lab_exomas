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


# Función para alinear
align_sample() {
    local DIR="$1"
    echo "$1"
    local SAMPLE=$(basename "$DIR")
    local FWD_FILE=$(find "$DIR" -type f \( -name "${SAMPLE}_R1.fq.gz" -o -name "${SAMPLE}_R1.fastq.gz" \))
    local REV_FILE=$(find "$DIR" -type f \( -name "${SAMPLE}_R2.fq.gz" -o -name "${SAMPLE}_R2.fastq.gz" \))
    local OUTPUT_DIR="$DIR"

    if [[ -f $FWD_FILE && -f $REV_FILE ]]; then
        echo "Procesando muestra: $SAMPLE"

        # Archivos intermedios
        local SAM_FILE="$OUTPUT_DIR/${SAMPLE}.sam"
        local BAM_FILE="$OUTPUT_DIR/${SAMPLE}.bam"
        local SORTED_BAM_FILE="$OUTPUT_DIR/${SAMPLE}_sorted.bam"
        local RG_BAM_FILE="$OUTPUT_DIR/${SAMPLE}_RG.bam"
        local MARKED_BAM_FILE="$OUTPUT_DIR/${SAMPLE}_deduped.bam"
        local METRICS_FILE="$OUTPUT_DIR/${SAMPLE}_metrics.txt"
        local FINAL_BAM_FILE_SORTED="$OUTPUT_DIR/${SAMPLE}_final_sorted.bam"

        # Montar la subcarpeta específica en Docker
        GATK_PATH="$GATK_PATH_BASE -v $DIR:/samples broadinstitute/gatk gatk"

        if [[ ! -f $FINAL_BAM_FILE_SORTED || ! -s $FINAL_BAM_FILE_SORTED ]]; then

            # Verificar si el archivo SAM ya existe y no está vacío
            if [[ -f $SAM_FILE && -s $SAM_FILE ]]; then
                echo "Archivo SAM para $SAMPLE ya existe y no está vacío. Omitiendo alineamiento."
            else
                # Verificar las rutas antes de ejecutar
                echo "Alineando $SAMPLE..."
                echo "Ejecutando: bwa mem -HMpP -v 3 -t 2 $REFERENCE $FWD_FILE $REV_FILE > $SAM_FILE"
                if [[ -f $REFERENCE && -f $FWD_FILE && -f $REV_FILE ]]; then
                    bwa mem -HMpP -v 3 -t 2 "$REFERENCE" "$FWD_FILE" "$REV_FILE" > "$SAM_FILE"
                else
                    echo "ERROR: Una o más rutas no son válidas."
                    echo "Referencia: $REFERENCE"
                    echo "Archivo Forward: $FWD_FILE"
                    echo "Archivo Reverse: $REV_FILE"
                    exit 1
                fi

                if [[ ! -s $SAM_FILE ]]; then
                    echo "ERROR: El archivo SAM generado para $SAMPLE está vacío. Revisa los datos de entrada o la referencia."
                    rm $SAM_FILE
                    exit 1
                fi
            fi

            # Verificar si el archivo BAM ya existe y no está vacío
            if [[ -f $BAM_FILE && -s $BAM_FILE ]]; then
                echo "Archivo BAM para $SAMPLE, el: $BAM_FILE, ya existe y no está vacío. Omitiendo alineamiento."
            else
                if [[ -f "$SAM_FILE" ]]; then
                    echo "Ejecutando: samtools view -S -b -h -@ 2 $SAM_FILE -o $BAM_FILE"
                    samtools view -S -b -h -@ 2 $SAM_FILE -o $BAM_FILE
                    if [[ ! -s $BAM_FILE ]]; then
                        echo "ERROR: El archivo BAM generado para $SAMPLE está vacío. Revisa los datos de entrada o la referencia."
                        rm $BAM_FILE
                        exit 1
                    else
                        echo "El archivo BAM se ha generado."
                    fi
                else
                    echo "ERROR: El archivo SAM falla."
                    echo "SAM FILE: $SAM_FILE"
                    exit 1
                fi
            fi

            # Ordenar BAM
            if [[ -f $SORTED_BAM_FILE && -s $SORTED_BAM_FILE ]]; then
                echo "Archivo BAM sorteado para $SAMPLE ya existe y no está vacío. Omitiendo."
            elif [[ -f $BAM_FILE ]]; then
                echo "Ejecutando: sorteando"
                $GATK_PATH SortSam -CREATE_INDEX true -INPUT "/samples/$(basename $BAM_FILE)" -OUTPUT "/samples/$(basename $SORTED_BAM_FILE)" -SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT

                if [[ ! -s $SORTED_BAM_FILE ]]; then
                    echo "ERROR: El archivo BAM sorteado generado para $SAMPLE está vacío. Revisa los datos de entrada o la referencia."
                    rm $SORTED_BAM_FILE
                    exit 1
                else
                    echo "El archivo BAM sorteado se ha generado."
                fi
            else
                echo "ERROR: Hay un problema con la ejecución del sort con GATK: $BAM_FILE"
            fi

            # Agregar Read Group
            if [[ -f $RG_BAM_FILE && -s $RG_BAM_FILE ]]; then
                echo "Archivo BAM deduplicado para $SAMPLE ya existe y no está vacío. Omitiendo."
            elif [[ -f $SORTED_BAM_FILE ]]; then
                echo "Agregando ReadGroup..."
                $GATK_PATH AddOrReplaceReadGroups \
                    I="/samples/$(basename $SORTED_BAM_FILE)" \
                    O="/samples/$(basename $RG_BAM_FILE)" \
                    RGID=$SAMPLE \
                    RGLB=lib2 \
                    RGPL=illumina \
                    RGPU=unit1 \
                    RGSM=$SAMPLE
                if [[ ! -s $RG_BAM_FILE ]]; then
                    echo "ERROR: El archivo BAM RG para $SAMPLE está vacío. Revisa los datos de entrada o la referencia."
                    rm $RG_BAM_FILE
                    exit 1
                else
                    echo "El archivo BAM deduplicado con Read Group se ha generado."
                fi
            else
                echo "ERROR: Hay un problema con $RG_BAM_FILE"
            fi

            # Marcar duplicados
            if [[ -f $MARKED_BAM_FILE && -s $MARKED_BAM_FILE && -f $METRICS_FILE && -s $METRICS_FILE ]]; then
                echo "Archivo BAM marcado con duplicados para $SAMPLE ya existe y no está vacío. Omitiendo."
            elif [[ -f $RG_BAM_FILE ]]; then
                echo '''$GATK_PATH MarkDuplicates -CREATE_INDEX true -INPUT "/samples/$(basename $RG_BAM_FILE)" -OUTPUT "/samples/$(basename $MARKED_BAM_FILE)" -M "/samples/$(basename $METRICS_FILE)" -VALIDATION_STRINGENCY STRICT --REMOVE_DUPLICATES true'''
                $GATK_PATH MarkDuplicates -CREATE_INDEX true -INPUT "/samples/$(basename $RG_BAM_FILE)" -OUTPUT "/samples/$(basename $MARKED_BAM_FILE)" -M "/samples/$(basename $METRICS_FILE)" -VALIDATION_STRINGENCY STRICT --REMOVE_DUPLICATES true

                if [[ ! -s "$MARKED_BAM_FILE" ]]; then
                    echo "ERROR: El archivo BAM deduplicado generado para $SAMPLE está vacío. Revisa los datos de entrada o la referencia."
                    rm $MARKED_BAM_FILE
                    exit 1
                else
                    echo "El archivo BAM deduplicado se ha generado."
                fi
            else
                echo "ERROR: Hay un problema con $RG_BAM_FILE"
            fi

            # Ordenar BAM final
            if [[ -f $FINAL_BAM_FILE_SORTED && -s $FINAL_BAM_FILE_SORTED ]]; then
                echo "Archivo BAM final ordenado para $SAMPLE ya existe y no está vacío. Omitiendo."
            elif [[ -f $MARKED_BAM_FILE ]]; then
                echo "Ejecutando: samtools sort -@ 2 -o $FINAL_BAM_FILE_SORTED $MARKED_BAM_FILE"
                samtools sort -@ 2 -o "$FINAL_BAM_FILE_SORTED" "$MARKED_BAM_FILE"
                ## creando el índice
                samtools index "$FINAL_BAM_FILE_SORTED"
                if [[ ! -s $FINAL_BAM_FILE_SORTED ]]; then
                    echo "ERROR: El archivo BAM final ordenado generado para $SAMPLE está vacío. Revisa los datos de entrada o la referencia."
                    rm "$FINAL_BAM_FILE_SORTED"
                    exit 1
                else
                    echo "El archivo BAM final ordenado se ha generado."
                fi
            else
                echo "ERROR: Hay un problema con $MARKED_BAM_FILE"
            fi

            ## remover archivos temporales
            echo "Eliminando archivos temporales. Solo nos quedamos con el deduped"
            echo "Ejecutando: rm $SAM_FILE $BAM_FILE $SORTED_BAM_FILE $RG_BAM_FILE $MARKED_BAM_FILE $METRICS_FILE"
            rm -f $SAM_FILE $BAM_FILE $SORTED_BAM_FILE $RG_BAM_FILE $MARKED_BAM_FILE $METRICS_FILE

        else
            echo "Esta alineado completamente, vamos a ver si no existen archivos temporales ... "
            if [[ -f $SAM_FILE || -f $BAM_FILE || -f $SORTED_BAM_FILE ||  -f $RG_BAM_FILE || -f $MARKED_BAM_FILE || -f $METRICS_FILE ]]; then
                echo "Existen los archivos temporales, vamos a eliminarlos"
                rm -f $SAM_FILE $BAM_FILE $SORTED_BAM_FILE $RG_BAM_FILE $MARKED_BAM_FILE $METRICS_FILE
            else
                echo "Todo correcto"
            fi
        fi
    fi
}

export -f process_sample  
export REFERENCE
export THREADS
export MAIN_DIR
export GATK_PATH_BASE  


# Procesar cada carpeta en paralelo
# Caso 1: Procesar un directorio específico
if [[ -n "$DIR" ]]; then
    if [[ -d "$DIR" ]]; then
        echo "Procesando directorio específico: $DIR"
        perform_fastqc "$MAIN_DIR"
        align_sample "$MAIN_DIR"
        
    else
        echo "El directorio $DIR no existe."
        exit 1
    fi
fi

# Caso 2: Buscar y procesar muestras específicas
if [[ -n "$SAMPLE_NAME" ]]; then
    for SAMPLE_DIR in $(find "$MAIN_DIR" -type d -name "$SAMPLE_NAME"); do
        echo "Procesando muestra específica: $SAMPLE_DIR"
        perform_fastqc "$SAMPLE_DIR"
        align_sample "$SAMPLE_DIR"
    done
fi
