from Analysis import Analysis

An = Analysis('assembly.gsd', 'assembly.map')


frame_set = [list(range(i * 100 , i * 100 + 3,1)) for i in range(1,3,1)]


frames = [f for frames in frame_set for f in frames]

for i in frames:
	print(i, An.average_connections(i, bound=False), An.invading_mers(i, bound=False), An.average_connections(i, bound=True), An.invading_mers(i, bound=True))
print(An.average_hybridization_lifetime(frames))
print(An.average_hybridization_lifetime(frames, bound=False ))


An.graph_rdf_combo_overlay(frame_set, cut_off=21)
