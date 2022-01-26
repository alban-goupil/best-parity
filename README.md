# Recherche des meilleures parités

## Objectifs

Étant donné un code de Gray d'une modulation d'ordre
supérieure comme des M-QAM ou des M-ASK, ils existent des
parités pour une longueur de code correcteur donnée qui
optimisent le spectre des distances euclidiennes. Leur
recherche est cependant difficile calculatoirement car le
domaine de recherche reste important. Les programmes de ce
répertoire utilisent des symétries et des coupures dans les
arbres de recherche pour limiter l'espace de recherche.

## Les programmes

Le répertoire [src/](./src/) contient les programmes
principaux avec une documentation dans les fichiers C.

Le répertoire [bin/](./src/) contient des exemples de
scripts qui lancent une recherche avec la possibilité de
paralléliser les processus avec le logiciel GNU parallel.

