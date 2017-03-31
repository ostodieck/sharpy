import sys
import sharpy.utils.cout_utils as cout


def read_settings():
    if len(sys.argv) == 1:
        cout.cout_wrap('Running SHARPy using the default settings file:')
        case_settings = './test.solver.txt'
        cout.cout_wrap('%s is the main settings file' % case_settings)
    elif len(sys.argv) == 2:
        case_settings = sys.argv[1]
        cout.cout_wrap('Running SHARPy using the settings file: %s' % case_settings)
    else:
        cout.cout_wrap('*** Too many arguments, only the first one will be used')
        case_settings = sys.argv[1]
        cout.cout_wrap('Running SHARPy using the settings file: %s' % case_settings)

    settings = parse_settings(case_settings)
    return settings


def parse_settings(file):
    import os
    from presharpy.utils.settings import load_config_file
    settings = load_config_file(os.path.realpath(file))
    settings['SHARPy']['flow'] = settings['SHARPy']['flow'].split('\n')
    return settings