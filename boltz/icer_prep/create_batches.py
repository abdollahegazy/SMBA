import math
from pathlib import Path

ligands = list(Path("../input").rglob("run.sh"))
batch_size = 1000 #max ICER jobs i believe
print(len(ligands))
print(math.ceil(len(ligands)/batch_size))

for i in range(math.ceil(len(ligands)/batch_size)):
    with open(f"batches/batch_{i+1}.txt","w") as f:
        for lig in ligands[i*batch_size:(i+1)*batch_size]:
            f.write(str(lig.parent) + "\n")
