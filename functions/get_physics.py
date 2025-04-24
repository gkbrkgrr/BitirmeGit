import pandas as pd

sets_mapping = pd.DataFrame([
    {"set_name": "SET1", "physics": "mp4_pbl1"},
    {"set_name": "SET2", "physics": "mp4_pbl5"},
    {"set_name": "SET3", "physics": "mp4_pbl7"},
    {"set_name": "SET4", "physics": "mp6_pbl1"},
    {"set_name": "SET5", "physics": "mp6_pbl5"},
    {"set_name": "SET6", "physics": "mp6_pbl7"},
    {"set_name": "SET7", "physics": "mp38_pbl1"},
    {"set_name": "SET8", "physics": "mp38_pbl5"},
    {"set_name": "SET9", "physics": "mp38_pbl7"},
    {"set_name": "SET10", "physics": "mp10_pbl1"},
    {"set_name": "SET11", "physics": "mp10_pbl5"},
    {"set_name": "SET12", "physics": "mp10_pbl7"},
    {"set_name": "SET13", "physics": "mp24_pbl1"},
    {"set_name": "SET14", "physics": "mp24_pbl5"},
    {"set_name": "SET15", "physics": "mp24_pbl7"},
])

def get_mp_and_pbl(set_name):
    physics = sets_mapping.loc[sets_mapping["set_name"] == set_name, "physics"].iloc[0]
    mp =  physics.split("_")[0].strip("mp")
    pbl = physics.split("_")[1].strip("pbl")
    return mp, pbl

