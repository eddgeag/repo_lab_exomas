#!/bin/bash

# Verifica que se haya proporcionado una ruta
if [ -z "$1" ]; then
  echo "Por favor, proporciona un path relativo donde están los archivos BAM."
  exit 1
fi

# Define la ruta proporcionada
relative_path="$1"

# Itera sobre cada archivo BAM en la ruta proporcionada y todas sus subcarpetas
find "$relative_path" -type f -name "*.bam" | while read -r bam_file; do
  # Verifica si el archivo BAM existe (la verificación ya se hace por `find`)
  if [ ! -e "$bam_file" ]; then
    echo "No se encontraron archivos BAM en el directorio especificado."
    exit 1
  fi

  # Ejecuta samtools quickcheck en el archivo BAM
  if ! samtools quickcheck "$bam_file"; then
    echo "Error en el archivo BAM: $bam_file"
  fi

  # Obtiene el nombre del archivo sin la extensión
  base_name="${bam_file%.bam}"

  # Elimina los archivos .bai que no correspondan al archivo BAM en cuestión y que no tengan el formato deseado
  find "$relative_path" -type f -name "*.bai" | while read -r bai_file; do
    bai_base_name="${bai_file%.bai}"

    # Comprobación para evitar eliminar los archivos con el formato DX0YZ-NM_final_sorted.bam.bai
    if [[ "$bai_file" =~ DX0[0-9]{2}-[0-9]{2}_final_sorted\.bam\.bai ]]; then
      continue
    fi

    # Elimina los archivos BAI que no coincidan con el archivo BAM actual
    if [[ "$bai_base_name" != "$base_name" ]]; then
      echo "Eliminando archivo BAI que no corresponde: $bai_file"
      rm  -f "$bai_file"
    fi
  done
done
