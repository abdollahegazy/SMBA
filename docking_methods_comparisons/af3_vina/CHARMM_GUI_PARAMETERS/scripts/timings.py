import os
import glob

COOMPLEXES_PATH = "../complexes"

print("Simulation Timings Check")
print("=" * 50)
print()

# Find all directories at level 3 (species/protein/ligand)
for root, dirs, files in os.walk(COOMPLEXES_PATH):
    depth = root.replace(COOMPLEXES_PATH, '').count(os.sep)
    
    # Only process directories at depth 3 (ligand level)
    if depth == 3:
        # Check both AF3_sim and MD_sim
        for sim_type in ["AF3_sim", "MD_sim"]:
            sim_dir = os.path.join(root, sim_type)
            
            if not os.path.isdir(sim_dir):
                continue
            
            # Skip toppar directories
            if os.path.basename(sim_dir) == "toppar":
                continue
            
            # Find all .out files
            out_files = glob.glob(os.path.join(sim_dir, "*.out"))
            
            if not out_files:
                print(f"❌ {sim_dir}")
                print(f"   No .out file found")
                print()
                continue
            
            # Get the latest .out file by modification time
            latest_out = max(out_files, key=os.path.getmtime)
            
            # Look for the WRITING VELOCITIES line
            step = None
            try:
                with open(latest_out, 'r') as f:
                    for line in f:
                        if "WRITING VELOCITIES TO OUTPUT FILE AT STEP" in line:
                            step = int(line.strip().split()[-1])
            except Exception as e:
                print(f"⚠️  {sim_dir}")
                print(f"   Error reading {os.path.basename(latest_out)}: {e}")
                print()
                continue
            
            if step is None:
                print(f"⚠️  {sim_dir}")
                print(f"   No velocity write found in {os.path.basename(latest_out)}")
                print()
                continue
            
            # Calculate nanoseconds: step * 2e-6
            ns = step * 2e-6
            
            print(f"✅ {sim_dir}")
            print(f"   Step: {step:,} → {ns:.2f} ns")
            print()