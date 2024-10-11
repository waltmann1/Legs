from __future__ import division
import numpy as np
import gsd.hoomd
import gsd.fl
import numpy as np
import numpy.linalg as la
import copy as cp
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 22})
import hoomd.data as hbox
from mpl_toolkits.mplot3d import Axes3D
import math as m

class Analysis2(object):

    def __init__(self, gsd_name, map_name, cut=.6, use_template=True, radius_approx=None):
        f = gsd.fl.open(name=gsd_name, mode='rb', application='', schema='hoomd',
                        schema_version=[1, 0])
        self.trajectory = gsd.hoomd.HOOMDTrajectory(f)
        self.invaders = None
        self.binders = None
        self.all_binders = None
        self.all_invaders = None
        self.arms = None
        self.xs = None
        self.ps = None
        self.temps = None
        self.read_map(map_name)
        self.frames = []
        self.topology = []
        self.connections = []
        self.bound = []
        self.box = self.trajectory.read_frame(0).configuration.box[:3]
        self.box_o = hbox.boxdim(Lx=self.box[0], Ly=self.box[1], Lz=self.box[2])
        self.cut =cut
        self.use_template = use_template
        self.radius_approx = radius_approx

    def read_map(self, map_name):

        f = open(map_name)
        data = f.readlines()
        self.binders = []
        self.all_binders = []
        self.invaders = []
        self.all_invaders = []
        self.arms = []
        self.xs = []
        self.ps = []
        self.temps = []
        on = 1000
        symbols = ['i', 'b', 'p', 'x', 'a', 't']
        mer_binders = []
        mer_invaders = []
        mer_arms = []
        mer_xs = []
        mer_ps = []
        mer_temps = []
        for line in data:
            s = line.split()
            for bit in s:
                if bit.isalpha():
                    on = symbols.index(bit)
                elif symbols[on] == 'i':
                    mer_invaders.append(int(bit))
                    self.all_invaders.append(int(bit))
                elif symbols[on] == 'b':
                    mer_binders.append(int(bit))
                    self.all_binders.append(int(bit))
                elif symbols[on] == 'p':
                    mer_ps.append(int(bit))
                elif symbols[on] == 'x':
                    mer_xs.append(int(bit))
                elif symbols[on] == 'a':
                    mer_arms.append(int(bit))
                elif symbols[on] == 't':
                    mer_temps.append(int(bit))
            if symbols[on] != 't':
                self.invaders.append(mer_invaders)
                mer_invaders = []
                self.binders.append(mer_binders)
                mer_binders = []
                self.ps.append(mer_ps)
                mer_ps = []
                self.xs.append(mer_xs)
                mer_xs = []
                self.arms.append(mer_arms)
                mer_arms = []
            else:
                self.temps.append(mer_temps)
                mer_temps = []

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
                            if la.norm(np.subtract(c_pos, np.add(frame.particles.position[d], np.multiply(self.box, frame.particles.image[d])))) < self.cut:
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
            #realpos = self.box_o.min_image(tuple(np.add(part, np.multiply(self.box, frame.particles.image[ind]))))
            realpos = frame.particles.image[ind]
            if ind > last:
                chain_pos.append(realpos)

        return chain_pos

    def template_removed(self, frame):

        if self.temps[0][0] == len(frame.particles.position) or not self.use_template:
            return True
        return False

    def bound_mers(self, frame_number, radius_approx=None):
        frame = self.trajectory.read_frame(frame_number)

        if self.template_removed(frame):
            return [0 for _ in range(len(self.invaders))]

        bound = []
        for i in range(len(self.invaders)):
            b = -1
            #position = self.box_o.min_image(tuple(np.add(frame.particles.position[i],
            #                                             np.multiply(self.box, frame.particles.image[i]))))
            position = frame.particles.position[i]
            for mer_ind, array in enumerate(self.temps):
                if radius_approx is not None:
                    pos = np.average(self.poly_positions([frame.particles.position[x] for x in array]), axis=0)
                    #print(pos)
                    if self.distance_image(position, pos) < 5 + radius_approx:
                        b = mer_ind
                        #print(mer_ind, array[0], frame_number, self.distance_image(position, pos), 5 + radius_approx)
                else:
                    for part_ind in array:
                        #pos = self.box_o.min_image(tuple(np.add(frame.particles.position[part_ind],
                        #             np.multiply(self.box, frame.particles.image[part_ind]))))
                        pos = frame.particles.position[part_ind]
                        #if la.norm(np.subtract(pos, position)) < 5:
                        if self.distance_image(position, pos) < 5:
                            b = mer_ind
            bound.append(b)
        #quit()
        return bound

    def read_frame(self, frame):

        if frame in self.frames:
            return
        else:
            self.frames.append(frame)
            top, connect = self.mer_topology(frame)
            self.topology.append(top)
            self.connections.append(connect)
            self.bound.append(self.bound_mers(frame, radius_approx=self.radius_approx))

    def distance_image(self, one, two):

        x_dist = np.abs(one[0] - two[0])
        if x_dist > self.box[0]/2:
            x_dist = self.box[0] - x_dist

        y_dist = np.abs(one[1] - two[1])
        if y_dist > self.box[1]/2:
            y_dist = self.box[1] - y_dist
        z_dist = np.abs(one[2] - two[2])
        if z_dist > self.box[2]/2:
            z_dist = self.box[2] - z_dist

        return la.norm([x_dist, y_dist, z_dist])

    def undo_pbc(self, one, two):

        dist = one[0] - two[0]
        #if dist > self.box[0]/2:
        if dist > self.box[0]/2:
            two[0] += self.box[0]
        elif dist < - self.box[0]/2:
            two[0] -= self.box[0]

        dist = one[1] - two[1]
        if dist > self.box[1]/2:
            two[1] += self.box[1]
        elif dist < -self.box[1]/2:
            two[1] -= self.box[1]

        dist = one[2] - two[2]
        if dist > self.box[2]/2:
            two[2] += self.box[2]
        elif dist < -self.box[2]/2:
            two[2] -= self.box[2]
        return two

    def polymer_end_to_end_distance(self, frame):

        snap = self.trajectory.read_frame(frame)
        index = 0
        while snap.particles.typeid[index] != snap.particles.types.index('qPm'):
            index += 1
        start = snap.particles.position[index]
        end = snap.particles.position[len(snap.particles.position) - 1]
        return self.distance_image(start, end)

    def poly_positions(self, pos):

        for i in range(len(pos) -1):
            pos[i+1] = self.undo_pbc(pos[i], pos[i+1])
        return pos


    def polymer_all_dists_squared(self, frame):
        snap = self.trajectory.read_frame(frame)
        index = 0
        while snap.particles.typeid[index] != snap.particles.types.index('qPm'):
            index += 1

        positions = snap.particles.position[index:]
        positions = self.poly_positions(positions)
        distances = []
        average_position = np.mean(positions, axis=0)
        for i in range(len(positions)):
            distances.append(self.distance_image(average_position, positions[i]))

        return [dist **2 for dist in distances]


    def polymer_radial_distribution(self, frame):
        snap = self.trajectory.read_frame(frame)
        index = 0
        cp_index = 0
        while snap.particles.typeid[index] != snap.particles.types.index('qPm'):
            index += 1
            if snap.particles.body[index] != -1:
                cp_index = index

        positions = list(snap.particles.position[:cp_index + 1])
        for new_pos in snap.particles.position[index:]:
            positions.append(new_pos)
        positions = self.poly_positions(positions)
        distances = []

        capsid_center = np.mean(positions[:cp_index + 1], axis=0)
        for pos in positions[cp_index+1:]:
            distances.append(self.distance_image(capsid_center, pos))

        return distances


    def average_connections(self, frame, bound=True):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        top = self.topology[frame_idx]
        b = self.bound[frame_idx]

        if bound:
            return np.average([len(cons) for ind,cons in enumerate(top) if len(cons) > 0 and b[ind] > -1])
        else:
            return np.average([len(cons) for ind,cons in enumerate(top) if len(cons) > 0 and b[ind] == -1])

    def invading_mers(self, frame, bound=True):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        top = self.topology[frame_idx]
        b = self.bound[frame_idx]

        if bound:
            return np.sum([1 for ind, cons in enumerate(top) if len(cons) > 0 and b[ind] > -1])
        else:
            return np.sum([1 for ind, cons in enumerate(top) if len(cons) > 0 and b[ind] == -1])

    def bound_correct(self, value, bound):

        if not bound:
            return value == -1
        else:
            return value > -1

    def bound_6nn_list(self, frame):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        b = self.bound[frame_idx]

        snap = self.trajectory.read_frame(frame)

        on_zero = [ind for ind, bound in enumerate(b) if bound == 0]

        list = [[] for _ in range(len(on_zero))]
        dist_list = [[] for _ in range(len(on_zero))]

        for oz_ind, mer_ind in enumerate(on_zero):
            for i in range(0, len(on_zero)):
                if i != oz_ind:
                    mer_ind2 = on_zero[i]
                    self.add_to_list_sorted(list[oz_ind], dist_list[oz_ind], mer_ind, mer_ind2, snap)
        return list,dist_list, on_zero

    def coordination_list(self, frame, cut=11):

        list,dist_list, on_zero = self.bound_6nn_list(frame)

        fos = []
        snap = self.trajectory.read_frame(frame)

        for ind, sub_list in enumerate(dist_list):

            ave = np.average(sub_list[:5])
            std = np.std(sub_list[:5])
            #print(sub_list, ind)
            if sub_list[5] < cut:
                fos.append(6)
            else:
                fos.append(5)

        positions = [self.min_image_position(mer_ind, snap) for mer_ind in on_zero]


        return positions, fos




    def min_image_position(self, mer_ind, snap):

        #return self.box_o.min_image(tuple(np.add(snap.particles.position[mer_ind],
                                #                            np.multiply(snap.particles.image[mer_ind], self.box))))

        return snap.particles.position[mer_ind]

    def add_to_list_sorted(self, mslist,msdlist, mer_ind, mer_ind2, snap):

        pos_two = self.min_image_position(mer_ind2, snap)
        pos_one = self.min_image_position(mer_ind, snap)
        #dist = la.norm(np.subtract(pos_one, pos_two))
        dist = self.distance_image(pos_one, pos_two)
        #print(len(msdlist), len(mslist))
        if len(mslist) == 0:
            mslist.append(mer_ind2)
            msdlist.append(dist)
            #print(len(msdlist), len(mslist))
            #print("")
            return
        else:
            for nn_rank, mer in enumerate(mslist):
                dist2 = msdlist[nn_rank]
                if dist < dist2:
                    mslist.insert(nn_rank, mer_ind2)
                    msdlist.insert(nn_rank, dist)
                    if len(mslist) == 7:
                        #print(mslist)
                        #print(msdlist)
                        mslist.remove(mslist[6])
                        msdlist.remove(msdlist[6])
                    #print(len(msdlist), len(mslist))
                    #print("")
                    return
            if len(mslist) < 6:
                mslist.append(mer_ind2)
                msdlist.append(dist)
                #print(len(msdlist), len(mslist))
                #print("")
                return



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
                    if self.bound_correct(bound_array[mer_ind], bound):
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
                if self.bound_correct(bound_array[mer_ind], bound):
                    for connect in x:
                        if connect in active_connections:
                            length[all_connections.index(connect)] += 1
                        else:
                            all_connections.append(connect)
                            length.append(1)
                            active_connections.append(connect)

        #print(length)
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

        mer_pos = [np.add(snap.particles.position[i], np.multiply(self.box, snap.particles.image[i])) for i in range(len(self.binders))]
        #mer_pos = [self.box_o.min_image(tuple(that)) for that in mer_pos]
        dists = []

        for ind, pos in enumerate(mer_pos):
            for mer_ind in range(ind + 1, len(mer_pos)):
                #dist = la.norm(np.subtract(pos, mer_pos[mer_ind]))
                dist = self.distance_image(pos, mer_pos[mer_ind])
                if bound[ind] == bound[mer_ind] and bound[ind] != -1 and dist < max(self.box)/2:
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
            #ax1.plot(bin_middles, rdf_hist)
        plt.legend()
        plt.show()

    def graph_rdf_combo(self, pdf_hist, cut_off=1000):

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

    def number_bound_distribution(self, frame):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        bound = self.bound[frame_idx]

        result = []

        for i in range(8):
            result.append(bound.count(i))

        return result


    def graph_polymer_rmsd(self, set_frames, show=True):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        #ax1.set_title("Frame " + str(frames))

        ax1.set_xlabel('time')
        ax1.set_ylabel('rmsd')

        times = [frame[0] for frame in set_frames]
        final_average = []
        final_std = []
        for frames in set_frames:
            all_dists = []
            for frame in frames:
                all_dists.append(self.polymer_end_to_end_distance(frame)**2)
            final_average.append(np.sqrt(np.mean(all_dists)))
            final_std.append(np.sqrt(np.std(all_dists)))


        ax1.errorbar(times, final_average, yerr=final_std)
        plt.legend()
        if show:
            plt.show()

    def graph_polymer_rg(self, set_frames, save_name=None, write=False, show=True):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        #ax1.set_title("Frame " + str(frames))

        ax1.set_xlabel('time')
        ax1.set_ylabel('Rg')

        times = [frame[0] for frame in set_frames]
        final_average = []
        for frames in set_frames:
            all_dists = []
            for frame in frames:
                all_dists.append(self.polymer_all_dists_squared(frame))
            all_dists = sum(all_dists, [])
            final_average.append(np.sqrt(np.sum(all_dists)/(len(all_dists))))

        ax1.plot(times, final_average)
        plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.pdf'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
            if write:
                save_name = save_name[:-4] + ".txt"
                f = open(save_name, 'w')
                for i in range(len(times)):
                    f.write(str(times[i]) + "  " + str(final_average[i]) + '\n')
        if show:
            plt.show()

    def graph_polymer_rdf_combo_overlay(self, set_frames, cut_off=1000, save_name=None, write_hist=False, show=True):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        #ax1.set_title("Frame " + str(frames))

        ax1.set_xlabel('r')
        ax1.set_ylabel('')

        if save_name is not None and write_hist:
            f = open(save_name + '.txt', 'w')

        for frames in set_frames:
            all_dists = []
            flat = []
            for frame in frames:
                all_dists.append(self.polymer_radial_distribution(frame))
            for dists in all_dists:
                for dist in dists:
                    flat.append(dist)

            bins= int((max(flat) - min(flat))*10)
            rdf_hist, rbe = np.histogram(flat, bins=bins)
            #print("len of flat", frame, len(flat))
            bin_middles = [(rbe[i] + rbe[i + 1]) / 2 for i in range(len(rbe) - 1)]

            #rdf_hist = [rdf_hist[i] / (4 * np.pi * bin_middles[i] * bin_middles[i]) for i in range(len(rdf_hist))]

            total = np.sum(rdf_hist)
            rdf_hist = [p / len(frames) for p in rdf_hist]

            ax1.plot(bin_middles, rdf_hist, label="Frame "+str(frames[0]))
            if save_name is not None and write_hist:
                f.write("Frame " + str(frames[0]) + '\n')
                for i in range(len(bin_middles)):
                    f.write(str(bin_middles[i]) + " " + str(rdf_hist[i]) + "\n")
        plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.pdf'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        if show:
            plt.show()

    def graph_rdf_combo_overlay(self, set_frames, cut_off=1000, save_name=None, write_hist=False, show=True):

        fig = plt.figure()

        ax1 = fig.add_subplot(111)

        ax1.set_title("Pair Distribution of VP1 Centers")

        ax1.set_xlabel('r (nm)')
        ax1.set_ylabel('Probability')

        if save_name is not None and write_hist:
            f = open(save_name + '.txt', 'w')
        count = 0
        for frames in set_frames:
            count += 1
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

            #ax1.plot(bin_middles, rdf_hist, label="Frame "+str(frames[0]))
            ax1.plot(bin_middles, np.divide(rdf_hist, np.sum(rdf_hist)),
                     label="t = " + str(count) + "/" + str(len(set_frames)),alpha=.3 + .7 * (count==2),
                     linewidth=1 + (count==2))
            #ax1.plot(bin_middles, np.divide(rdf_hist, np.sum(rdf_hist)))
            #ax1.plot(bin_middles, np.divide(rdf_hist, np.sum(rdf_hist)), label="9nm template")
            if save_name is not None and write_hist:
                f.write("Frame "+str(frames[0]) + '\n')
                for i in range(len(bin_middles)):
                    f.write(str(bin_middles[i]) + " " + str(rdf_hist[i]) + "\n")
        plt.legend(fontsize=10)
        ax1.set_ylim([0, .055])
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        if show:
            plt.show()


    def relative_spherical_positions(self, frames):

        centers = []
        total_positions = []
        for frame in frames:
            positions = []
            frame = self.trajectory.read_frame(frame)
            count = 1
            positions.append(frame.particles.position[0])
            while frame.particles.body[count] > frame.particles.body[count-1]:
                positions.append(frame.particles.position[count])
                count += 1
            centers.append(positions[-1])
            positions = positions[:-1]
            total_positions.append(positions)
        return self.spherical_average(total_positions, centers)

    def spherical_average(self, positions, centers):

        print("centers shape", np.array(centers).shape)
        print("positions shape", np.array(positions).shape)
        total = np.zeros_like(np.array(positions[0]))

        for ind, frame in enumerate(positions):
            for pos_ind, pos in enumerate(frame):
                pos = self.undo_pbc(centers[ind], pos)
                pos = np.subtract(pos, centers[ind])
                sph_pos = self.cart2sphA(pos)
                total[pos_ind] = np.add(total[pos_ind], sph_pos)

        average = np.divide(total, len(centers))

        print("total shape", total.shape)
        print("average shape", average.shape)
        print("average", average)

        average_cart_pos = [self.sph2cartA(pos[0], pos[1], pos[2]) for pos in average]
        print("average", average_cart_pos)
        return average_cart_pos

    def cart2sph(self, x, y, z):
        XsqPlusYsq = x ** 2 + y ** 2
        r = m.sqrt(XsqPlusYsq + z ** 2)  # r
        #phi = m.atan2(z, m.sqrt(XsqPlusYsq))  # phi
        phi = m.acos(z/r)
        theta = m.atan2(y, x)  # theta
        return r, theta, phi

    def cart2sphA(self, pt):
        r, theta, phi = self.cart2sph(pt[0], pt[1], pt[2])
        return np.array([r, theta, phi])

    def sph2cartA(self, r, theta, phi):

        x = r * np.cos(theta) * np.sin(phi)
        y = r * np.sin(theta) * np.sin(phi)
        if phi < 0:
            phi += np.pi
        z = r * np.cos(phi)
        return [x, y, z]


    def graph_averaged_coord_sphere(self, frames, save_name=None):

        on_zero = self.relative_spherical_positions(frames)
        list = [[] for _ in range(len(on_zero))]
        dist_list = [[] for _ in range(len(on_zero))]

        for oz_ind, mer_pos in enumerate(on_zero):
            for i in range(0, len(on_zero)):
                if i != oz_ind:
                    mer_pos2 = on_zero[i]
                    dist_list[i].append(self.distance_image(mer_pos, mer_pos2))
        fos = []
        for i in range(len(list)):
            dist_list[i].sort()
            if dist_list[i][4]> 13:
                fos.append(4)
            elif dist_list[i][5] > 13:
                fos.append(5)
            else:
                fos.append(6)


        xs = [pos[0] for pos in on_zero]
        ys = [pos[1] for pos in on_zero]
        zs = [pos[2] for pos in on_zero]
        colors = []
        labels = []

        for coord in fos:
            if coord == 6:
                colors.append('k')
                labels.append("6")
            elif coord == 5:
                colors.append('red')
                labels.append("5")
            else:
                colors.append("yellow")
                labels.append("4")

        fig = plt.figure()

        ax1 = fig.add_subplot(111, projection='3d')

        ax1.set_title("Black is 6 red is 5 yellow is 4")

        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')

        ax1.scatter(xs, ys, zs, c=colors, s=1000, label=labels)

        #plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.pdf'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()



    def graph_coord_sphere(self, frame, save_name=None, cut=11):

        fig = plt.figure()

        ax1 = fig.add_subplot(111, projection='3d')

        ax1.set_title("Black is 6 red is 5 yellow is other ")

        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        positions, fos = self.coordination_list(frame, cut=cut)
        xs = [pos[0] for pos in positions]
        ys = [pos[1] for pos in positions]
        zs = [pos[2] for pos in positions]
        colors = []
        labels = []

        for coord in fos:
            if coord == 6:
                colors.append('k')
                labels.append("6")
            elif coord == 5:
                colors.append('red')
                labels.append("5")
            else:
                colors.append("yellow")
                labels.append("other")

        ax1.scatter(xs, ys, zs, c=colors, s=1000, label=labels)

        #plt.legend()
        if save_name is not None:
            if save_name[-4] != ".":
                save_name += '.png'
            plt.savefig(save_name, bbox_inches='tight', pad_inches=.2)
        plt.show()

    def get_doubles(self, frame):

        if frame not in self.frames:
            self.read_frame(frame)
            frame_idx = len(self.frames) - 1
        else:
            frame_idx = self.frames.index(frame)

        connections = self.connections[frame_idx]

        doubles = 0

        for own_mer_ind, mer in enumerate(connections):
            for connect in mer:
                own_bind_ind = self.invaders[own_mer_ind].index(connect[0])
                opp_mer_ind, opp_invader = self.get_binder_invader_pair_index(connect[1])
                for back_connect in connections[opp_mer_ind]:
                    if opp_invader in back_connect:
                        own_binder = back_connect[1]
                        if own_binder == self.binders[own_mer_ind][own_bind_ind]:
                            doubles += 1
        return doubles

    def get_binder_invader_pair_index(self, binder_ind):

        for mer_ind, mer in enumerate(self.binders):
            for bind_ind, ind in enumerate(mer):
                if ind == binder_ind:
                    return mer_ind, self.invaders[mer_ind][bind_ind]
