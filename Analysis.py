from __future__ import division
import numpy as np
import gsd.hoomd
import gsd.fl
import numpy as np
import numpy.linalg as la
import copy as cp
from matplotlib import pyplot as plt

class Analysis(object):

    def __init__(self, gsd_name, map_name):
        f = gsd.fl.open(name=gsd_name, mode='rb', application='', schema='hoomd',
                        schema_version=[1, 0])
        self.trajectory = gsd.hoomd.HOOMDTrajectory(f)
        self.invaders = None
        self.binders = None
        self.all_binders = None
        self.all_invaders = None
        self.read_map(map_name)
        self.frames = []
        self.topology = []
        self.connections = []
        self.bound = []
        self.box = self.trajectory.read_frame(0).configuration.box[:3]

    def read_map(self, map_name):

        f = open(map_name)
        data = f.readlines()
        self.binders = []
        self.all_binders = []
        self.invaders = []
        self.all_invaders = []
        mer_binders = []
        mer_invaders = []
        for line in data:
            s = line.split()[1:]
            invader = True
            for bit in s:
                if bit.isalpha():
                    invader = False
                elif invader:
                    mer_invaders.append(int(bit))
                    self.all_invaders.append(int(bit))
                else:
                    mer_binders.append(int(bit))
                    self.all_binders.append(int(bit))
            self.invaders.append(mer_invaders)
            mer_invaders = []
            self.binders.append(mer_binders)
            mer_binders = []

    def mer_topology(self, frame):

        frame = self.trajectory.read_frame(frame)
        top = []
        cd_connections = []
        for ind, mer in enumerate(self.invaders):
            connections = []
            explicit_connections = []
            for c in mer:
                c_pos = np.add(frame.particles.position[c], np.multiply(self.box, frame.particles.image[c]))
                found = False
                for mer_ind, sites in enumerate(self.binders):
                    if not found:
                        for d in sites:
                            if la.norm(np.subtract(c_pos, np.add(frame.particles.position[d], np.multiply(self.box, frame.particles.image[d])))) < .6:
                                found = True
                                connections.append(mer_ind)
                                explicit_connections.append([c,d])

            top.append(connections)
            cd_connections.append(explicit_connections)

        return top, cd_connections

    def chain_positions(self, frame_number):

        last = self.invaders[-1][-1]
        chain_pos = []
        frame = self.trajectory.read_frame(frame_number)

        for ind, part in enumerate(frame.particles.position):
            realpos = np.add(part, np.multiply(self.box, frame.particles.image[ind]))
            if ind > last:
                chain_pos.append(realpos)

        return chain_pos

    def bound_mers(self, frame_number):
        frame = self.trajectory.read_frame(frame_number)

        bound = []
        chains = self.chain_positions(frame_number)
        for i in range(len(self.invaders)):
            b = 0
            position = np.add(frame.particles.position[i], np.multiply(self.box, frame.particles.image[i]))
            for pos in chains:
                if la.norm(np.subtract(pos, position)) < 5:
                    b = 1
            bound.append(b)
        return bound

    def read_frame(self, frame):

        if frame in self.frames:
            return
        else:
            self.frames.append(frame)
            top, connect = self.mer_topology(frame)
            self.topology.append(top)
            self.connections.append(connect)
            self.bound.append(self.bound_mers(frame))

    def average_connections(self, frame, bound=True):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        top = self.topology[frame_idx]
        b = self.bound[frame_idx]

        return np.average([len(cons) for ind,cons in enumerate(top) if len(cons) > 0 and b[ind] == bound])

    def invading_mers(self, frame, bound=True):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        top = self.topology[frame_idx]
        b = self.bound[frame_idx]

        return np.sum([1 for ind, cons in enumerate(top) if len(cons) > 0 and b[ind] == bound])

    def hybridization_lifetimes(self,frames, bound=True):
        """

        :param frames: assumed to be in order first to last
        :param bound: for the mers bound to the template or not bound
        :return:
        """

        hybrids = []

        for frame in frames:
            if frame not in self.frames:
                self.read_frame(frame)
                frame_idx = len(self.frames) - 1
            else:
                frame_idx = self.frames.index(frame)

            con = self.connections[frame_idx]
            hybrids.append(con)

        active_connections = []
        all_connections = []
        length = []

        for mer_ind, mer in enumerate(hybrids[0]):
            for connect in mer:
                    frame_idx = self.frames.index(frames[0])
                    bound_array = self.bound[frame_idx]
                    if bound_array[mer_ind] == bound:
                        all_connections.append(connect)
                        active_connections.append(connect)
                        length.append(1)


        for ind, f in enumerate(hybrids[1:]):
            now = [n_connect for mert in f for n_connect in mert]

            for con in active_connections:
                if con not in now:
                    active_connections.remove(con)

            for mer_ind, x in enumerate(f):
                frame_idx = self.frames.index(frames[ind + 1])
                bound_array = self.bound[frame_idx]
                if bound_array[mer_ind] == bound:
                    for connect in x:
                        if connect in active_connections:
                            length[all_connections.index(connect)] += 1
                        else:
                            all_connections.append(connect)
                            length.append(1)
                            active_connections.append(connect)

        print(length)
        return all_connections, length

    def average_hybridization_lifetime(self, frames, bound=True):

        all_connections, length = self.hybridization_lifetimes(frames, bound=bound)

        return(np.average(length))

    def bound_radial_distribution(self, frame, cut_off=1000):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        bound = self.bound[frame_idx]

        snap = self.trajectory.read_frame(frame)

        mer_pos = [np.add(snap.particles.position[i], np.multiply(self.box, snap.particles.image[i])) for i in range(len(bound)) if bound[i]]
        dists = []

        for ind, pos in enumerate(mer_pos):
            for mer_ind in range(ind + 1, len(mer_pos)):
                dist = la.norm(np.subtract(pos, mer_pos[mer_ind]))
                if dist < cut_off:
                    dists.append(dist)

        return dists

    def graph_rdf(self, frame, cut_off=1000):

        dists = self.bound_radial_distribution(frame, cut_off=cut_off)
        rdf_hist, rbe = np.histogram(dists, bins=100)
        bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

        pdf_hist = [rdf_hist[i] / (4 * np.pi * bin_middles[i] * bin_middles[i]) for i in range(len(rdf_hist))]

        total = np.sum(pdf_hist)

        pdf_hist = [p / total for p in pdf_hist]
        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frame))

        ax1.set_xlabel('r')
        ax1.set_ylabel('')
        ax1.plot(bin_middles, rdf_hist, label='RDF')
        plt.show()

    def graph_rdf_overlay(self, frames, cut_off=1000):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        #ax1.set_title("Frame " + str(frames))

        ax1.set_xlabel('r')
        ax1.set_ylabel('')
        for frame in frames:
            dists = self.bound_radial_distribution(frame, cut_off=cut_off)
            rdf_hist, rbe = np.histogram(dists, bins=100)
            bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

            #pdf_hist = [rdf_hist[i] / (4 * np.pi * bin_middles[i] * bin_middles[i]) for i in range(len(rdf_hist))]

            total = np.sum(rdf_hist)

            pdf_hist = [p / total for p in rdf_hist]

            ax1.plot(bin_middles, rdf_hist, label="Frame "+str(frame))
        plt.legend()
        plt.show()

    def graph_rdf_combo(self, frames, cut_off=1000):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Frame " + str(frames))

        ax1.set_xlabel('r')
        ax1.set_ylabel('')
        all_dists = []
        flat = []
        for frame in frames:
            all_dists.append(self.bound_radial_distribution(frame, cut_off=cut_off))
        for dists in all_dists:
                for dist in dists:
                        flat.append(dist)
        rdf_hist, rbe = np.histogram(flat, bins=100)
        bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

        #pdf_hist = [rdf_hist[i] / (4 * np.pi * bin_middles[i] * bin_middles[i]) for i in range(len(rdf_hist))]

        total = np.sum(rdf_hist)

        pdf_hist = [p / total for p in rdf_hist]

        ax1.plot(bin_middles, rdf_hist, label="Frame "+str(frame))
        plt.legend()
        plt.show()

    def graph_rdf_combo_overlay(self, set_frames, cut_off=1000):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        #ax1.set_title("Frame " + str(frames))

        ax1.set_xlabel('r')
        ax1.set_ylabel('')

        for frames in set_frames:
            all_dists = []
            flat = []
            for frame in frames:
                all_dists.append(self.bound_radial_distribution(frame, cut_off=cut_off))
            for dists in all_dists:
                for dist in dists:
                    flat.append(dist)

            rdf_hist, rbe = np.histogram(flat, bins=100)
            bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

            #pdf_hist = [rdf_hist[i] / (4 * np.pi * bin_middles[i] * bin_middles[i]) for i in range(len(rdf_hist))]

            total = np.sum(rdf_hist)
            rdf_hist = [p / len(frames) for p in rdf_hist]

            ax1.plot(bin_middles, rdf_hist, label="Frame "+str(frames[0]))
        plt.legend()
        plt.show()

