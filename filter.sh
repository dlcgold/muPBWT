#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Uso: $0 reference.bcf target.bcf"
  exit 1
fi

REF_BCF="$1"
TGT_BCF="$2"

OUT_REF="f_$1"
OUT_TGT="f_$2"

TMP_POS="positions.txt"
TMP_DISCORDANT="discordant_variants.bcf"

bcftools index "$1"
bcftools index "$2"

bcftools query -f '%CHROM\t%POS\n' "$TGT_BCF" | sort -u >"$TMP_POS"

bcftools view -T "$TMP_POS" "$REF_BCF" -Ob -o ref_temp.bcf
bcftools view -T "$TMP_POS" "$TGT_BCF" -Ob -o tgt_temp.bcf

bcftools index ref_temp.bcf
bcftools index tgt_temp.bcf

bcftools isec -c none -n=1 ref_temp.bcf tgt_temp.bcf -Ob -o "$TMP_DISCORDANT"

bcftools view -T ^"$TMP_DISCORDANT" ref_temp.bcf -Ob -o "$OUT_REF"
bcftools view -T ^"$TMP_DISCORDANT" tgt_temp.bcf -Ob -o "$OUT_TGT"

rm -f "$TMP_POS" "$TMP_DISCORDANT" ref_temp.bcf tgt_temp.bcf "$TMP_DISCORDANT".csi ref_temp.bcf.csi tgt_temp.bcf.csi

echo "Filtraggio completato! Output:"
echo " - Reference filtrato: $OUT_REF"
echo " - Target filtrato: $OUT_TGT"
