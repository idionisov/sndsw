import numpy as np
import pyhepmc
import uproot


def to_ntuple(
    input_file: str,
    output_file: str,
    tree: str = "nt"
):
    data = {
        "run": [], "event": [], "id": [], "generation": [],
        "E": [], "w": [],
        "x": [], "y": [], "z": [], "t": [], "px": [], "py": [], "pz": []
    }

    print(f"Reading {input_file}...")
    with pyhepmc.open(input_file) as f:
        for event in f:
            evt_num = event.event_number

            weight = event.weights[0] if event.weights else 1.0
            for particle in event.particles:
                vtx = particle.production_vertex
                if vtx:
                    pos = vtx.position
                    x, y, z, t = pos.x, pos.y, pos.z, pos.t
                else:
                    x, y, z, t = 0.0, 0.0, 0.0, 0.0

                mom = particle.momentum

                data["run"].append(0)
                data["event"].append(evt_num)
                data["id"].append(particle.pid)
                data["generation"].append(particle.status)
                data["w"].append(weight)

                data["E"].append(mom.e)
                data["px"].append(mom.px)
                data["py"].append(mom.py)
                data["pz"].append(mom.pz)

                data["x"].append(x)
                data["y"].append(y)
                data["z"].append(z)
                data["t"].append(t)

    print(f"Writing {len(data['x'])} particles to {output_file}...")

    numpy_data = {key: np.array(value) for key, value in data.items()}
    with uproot.recreate(output_file) as file:
        file[tree] = numpy_data

    print("Done!")


to_ntuple(input_file="/afs/cern.ch/user/i/idioniso/FORESEE/Models/DarkHiggs/model/events/test.hepmc", output_file="/eos/user/i/idioniso/hepmc.root")
