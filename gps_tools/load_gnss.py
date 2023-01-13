
import sys, os
from .file_io import config_io, io_nota, io_magnet_unr, io_usgs
from . import gps_input_vel_pipeline, vel_functions, gps_input_pipeline, utilities


def create_station_repo(root_config, refframe, proc_center, subnetwork=None):
    if proc_center == 'cwu':
        engine = cwu_proc_engine(root_config, refframe);
    elif proc_center == 'unr':
        engine = unr_proc_engine(root_config, refframe);
    elif proc_center == 'usgs':
        engine = usgs_proc_engine(root_config, refframe, subnetwork);
    elif proc_center == 'pbo':
        engine = pbo_proc_engine(root_config, refframe, proc_center);
    elif proc_center == 'nmt':
        engine = pbo_proc_engine(root_config, refframe, proc_center);
    else:
        print("Error! Invalid choice of network [pick one of pbo/cwu/unr/usgs]");
        sys.exit(0);
    return station_repo(engine)


class station_repo:
    # Loads up the general database functionality.
    # Contains one object called proc_engine, and a bunch of services.

    def __init__(self, proc_engine):
        self.proc_engine = proc_engine;
        print("Building database of stations using %s " % self.proc_engine.print_engine_metadata());

    def import_full_velfield(self):
        velfield = self.proc_engine.import_velfield();
        return velfield;

    def search_stations_by_circle(self, center, radius, basic_clean=False):
        """
        :param center: list, [lon, lat]
        :param radius: float, km
        :param basic_clean: bool.  If True, then keep long-lived stations and stations with low velocity uncertainties
        :returns: list of station_vels, list of floats [stations, distances]
        """
        velfield = self.proc_engine.import_velfield();
        if basic_clean:  # take long-lived stations and stations with low velocity uncertainties
            velfield = vel_functions.basic_clean_stations(velfield);
        close_stations, rad_distance = vel_functions.filter_to_circle(velfield, center, radius);
        return close_stations, rad_distance;

    def search_stations_by_box(self, bbox, basic_clean=False):
        """
        :param bbox: [lonW, lonE, latS, latN]
        :param basic_clean: bool.  If True, then keep long-lived stations and stations with low velocity uncertainties
        :returns: list of station_vels, list of floats [stations, distances]
        """
        velfield = self.proc_engine.import_velfield();
        if basic_clean:  # take long-lived stations and stations with low velocity uncertainties
            velfield = vel_functions.basic_clean_stations(velfield);
        close_stations = vel_functions.filter_to_bounding_box(velfield, bbox);
        return close_stations;

    def load_stations(self, station_name_list):
        """Load the full time series objects from a list of station names and a database"""
        data_obj_list, offset_obj_list, eq_obj_list = [], [], [];
        for station in station_name_list:
            [myData, offset_obj, eq_obj] = self.proc_engine.import_station(station);
            data_obj_list.append(myData);
            offset_obj_list.append(offset_obj);
            eq_obj_list.append(eq_obj);
        return [data_obj_list, offset_obj_list, eq_obj_list]


class cwu_proc_engine:
    # specific read/write functions. Each one has a matching load velfield, load station

    def __init__(self, root_config, refframe):
        self.root_config = root_config;   # string
        self.refframe = refframe;  # string
        self.file_params = config_io.read_config_file(self.root_config);

    def print_engine_metadata(self):
        return "Network: CWU, Refframe: %s " % self.refframe;

    def import_velfield(self):  # load the entire velocity field of CWU in a certain reference frame
        lookup_dict = gps_input_vel_pipeline.build_lookup_dictionary(self.file_params, 'cwu');  # read filenames
        myVelocities = io_nota.read_pbo_vel_file_format(lookup_dict["cwu_" + self.refframe]);
        return myVelocities;

    def import_station(self, station):
        filename = self.pre_screen_datasource_paths(station);
        [myData, offset_obj, eq_obj] = gps_input_pipeline.get_cwu(self.file_params, filename, station, 'cwu');
        return [myData, offset_obj, eq_obj];

    def pre_screen_datasource_paths(self, station):
        if self.refframe == "NA":
            filename = self.file_params['cwu']["directory"] + self.file_params['cwu']["gps_ts_dir"] + \
                       station + '.cwu.final_nam14.pos';
        elif self.refframe == "ITRF":  # ITRF
            filename = self.file_params['cwu']["directory"] + self.file_params['cwu']["gps_ts_dir"] + \
                       station + '.cwu.final_igs14.pos';
        else:
            print("Error! Reference frame doesn't match available ones [NA, ITRF]. Choose again"); sys.exit(1);
        utilities.check_if_file_exists(filename);
        return filename;


class unr_proc_engine:
    # specific read/write functions. Each one has a matching load velfield, load station

    def __init__(self, root_config, refframe):
        self.root_config = root_config;   # string
        self.refframe = refframe;  # string
        self.file_params = config_io.read_config_file(self.root_config);

    def print_engine_metadata(self):
        return "Network: UNR, Refframe: %s " % self.refframe;

    def import_velfield(self):
        lookup_dict = gps_input_vel_pipeline.build_lookup_dictionary(self.file_params, 'unr');  # read filenames
        myVelocities = io_magnet_unr.read_unr_vel_file(lookup_dict["unr_" + self.refframe],
                                                       self.file_params["unr"]["directory"] +
                                                       self.file_params["unr"]["coords_file"]);
        return myVelocities;

    def import_station(self, station_name):
        filename = self.pre_screen_datasource_paths(station_name);
        [myData, offset_obj, eq_obj] = gps_input_pipeline.get_unr(self.file_params, filename);  # UNR data format
        return [myData, offset_obj, eq_obj];

    def pre_screen_datasource_paths(self, station):
        if self.refframe == 'NA':
            filename = self.file_params["unr"]["directory"]+self.file_params["unr"]["gps_ts_dir"]+station+'.NA.tenv3';
        elif self.refframe == 'ITRF':
            filename = self.file_params["unr"]["directory"]+self.file_params["unr"]["gps_ts_dir"]+station+'.IGS14.tenv3'
        else:
            print("Error! Reference frame doesn't match available ones [NA, ITRF]. Choose again"); sys.exit(1);
        utilities.check_if_file_exists(filename);
        return filename;


class usgs_proc_engine:
    # specific read/write functions. Each one has a matching load velfield, load station

    def __init__(self, root_config, refframe, subnetwork):
        self.root_config = root_config;   # string
        self.refframe = refframe;  # string
        self.subnetwork = subnetwork;  # string
        self.file_params = config_io.read_config_file(self.root_config);

    def print_engine_metadata(self):
        return "Network: USGS, Subnetwork: %s, Refframe: %s " % (self.subnetwork, self.refframe);

    def import_velfield(self):
        lookup_dict = gps_input_vel_pipeline.build_lookup_dictionary(self.file_params, 'usgs', self.subnetwork);
        myVelocities = io_usgs.read_usgs_velfile(lookup_dict["usgs_" + self.refframe],
                                                 self.file_params["usgs"]["directory"] +
                                                 self.file_params["usgs"]["cache_file"]);
        return myVelocities;

    def import_station(self, station_name):
        filename, sub_network = self.pre_screen_datasource_paths(station_name);
        [myData, offset_obj, eq_obj] = gps_input_pipeline.get_usgs(self.file_params, filename, station_name,
                                                                   self.refframe, self.subnetwork);  # USGS data
        return [myData, offset_obj, eq_obj];

    def pre_screen_datasource_paths(self, station_name):
        # CHECK IF STATION EXISTS IN JUST ONE SUBNETWORK FOR CONVENIENCE
        if self.subnetwork == '':
            network_list = io_usgs.query_usgs_network_name(station_name, self.file_params['usgs']['directory'] +
                                                           self.file_params['usgs']['gps_ts_dir']);
            if len(network_list) == 1:
                self.subnetwork = network_list[0].split('/')[-1];
            else:
                print("ERROR! User must select one sub-network for USGS time series. Exiting. ");
                sys.exit(1);
        if self.refframe == "NA":
            filename = self.file_params["usgs"]["directory"]+self.file_params["usgs"]["gps_ts_dir"]+self.subnetwork + \
                       '/' + station_name.lower() + "_NAfixed.rneu";
        elif self.refframe == "ITRF":
            filename = self.file_params["usgs"]["directory"]+self.file_params["usgs"]["gps_ts_dir"]+self.subnetwork + \
                       '/' + station_name.lower() + "_ITRF2008.rneu";
        else:
            print("Error! Reference frame doesn't match available ones [NA, ITRF]. Choose again"); sys.exit(1);
        utilities.check_if_file_exists(filename);
        if os.path.isfile(filename):  # Determine if file is found on system. Provide helpful suggestions if not.
            print("Found file %s from datasource " % filename);
        else:  # If the file is not found on the system:
            print("Error!  Cannot find %s in database." % filename);
            print("The station %s could be found in the following sub_networks instead: " % station_name);
            io_usgs.query_usgs_network_name(station_name, self.file_params["usgs"]["directory"]);
            sys.exit(1);
        return filename, self.subnetwork;


class pbo_proc_engine:
    # specific read/write functions. Each one has a matching load velfield, load station

    def __init__(self, root_config, refframe, proccenter):
        self.root_config = root_config;   # string
        self.refframe = refframe;  # string
        self.file_params = config_io.read_config_file(self.root_config);
        self.proccenter = proccenter;

    def print_engine_metadata(self):
        return "Network: %s, Refframe: %s " % (self.proccenter, self.refframe);

    def import_velfield(self):  # load the entire velocity field of CWU in a certain reference frame
        lookup_dict = gps_input_vel_pipeline.build_lookup_dictionary(self.file_params, self.proccenter);
        myVelocities = io_nota.read_pbo_vel_file_format(lookup_dict[self.proccenter + "_" + self.refframe]);
        return myVelocities;

    def import_station(self, station):
        filename = self.pre_screen_datasource_paths(station);
        [myData, offset_obj, eq_obj] = gps_input_pipeline.get_pbo_type(self.file_params, filename, station,
                                                                       self.proccenter);
        return [myData, offset_obj, eq_obj];

    def pre_screen_datasource_paths(self, station):
        if self.refframe == "NA":
            filename = self.file_params[self.proccenter]["directory"]+self.file_params[self.proccenter]["gps_ts_dir"] \
                       + station + '.'+self.proccenter+'.final_nam08.pos';
        elif self.refframe == "ITRF":  # ITRF
            filename = self.file_params[self.proccenter]["directory"]+self.file_params[self.proccenter]["gps_ts_dir"] \
                       + station + '.'+self.proccenter+'.final_igs08.pos';
        else:
            print("Error! Reference frame doesn't match available ones [NA, ITRF]. Choose again"); sys.exit(1);
        utilities.check_if_file_exists(filename);
        return filename;
