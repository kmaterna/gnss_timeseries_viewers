"""
Functions for reading a variety of GNSS data formats, including:
  GNSS time series
  GNSS velocity fields
  GRACE and other hydrological models
"""

import collections, sys, os, configparser


def read_config_file(infile):
    """
    Reading system configuration into a dictionary of dictionaries.
    Upper level keys are for each database: [unr, pbo, cwu, etc.]

    :param infile: string, filename
    :returns: dictionary, which may contain other dictionaries
    """
    if not os.path.isfile(infile):
        print("Error! Data Config file %s not found on your machine. Must fix!" % infile);
        sys.exit(1);

    # Read the all important config file.
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(infile);
    config_section = config["py-config"];

    # Create a default dictionary so we can tolerate a config file with less-complete fields
    param_dict = collections.defaultdict(lambda: "Key Not Present In Config");   # dictionary of dictionaries
    for key in config_section.keys():
        param_dict[key] = config_section[key];
        if 'config' in key:
            database_name = key.split('_')[0];
            one_dictionary = read_one_database_config(config_section[key], 'data-config');
            one_dictionary['directory'] = "/".join(param_dict[key].split('/')[0:-1])+'/';  # the directory for data
            # one_dictionary['directory'] = os.path.join(list(param_dict[key].split('/')[0:-1]));  # WORK IN PROGRESS
            param_dict[database_name] = one_dictionary;
    return param_dict;


def read_one_database_config(configfile, sectionname):
    """
    :param configfile: string, filename
    :param sectionname: string
    :returns: dictionary, representing one database
    """
    config = configparser.ConfigParser()
    config.optionxform = str  # make the config file case-sensitive
    config.read(configfile);
    config_dictionary = dict(config[sectionname]);
    return config_dictionary;
