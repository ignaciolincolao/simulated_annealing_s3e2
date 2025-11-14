#!/usr/bin/env python
# encoding: utf-8

import os
from waflib.Configure import conf


def options(opt):
  opt.add_option('--paralelo', type='string', help='path to paralelo', dest='paralelo')

@conf
def check_paralelo(conf):
    # Ejemplo de comprobación de la biblioteca "paralelo"
    if conf.options.paralelo:
        includes_check = [conf.options.paralelo + '/include']
        libs_check = [conf.options.paralelo + '/lib']
    else:
        includes_check = ['/usr/local/include/paralelo']
        libs_check = ['/usr/local/lib']
    try:
        conf.start_msg('Checking for Paralelo includes')
        conf.find_file('SimulatedFactory.hpp', includes_check)
        conf.end_msg('ok')
        conf.start_msg('Checking for Paralelo libs')
        conf.find_file('libparalelo.so', libs_check)
        conf.end_msg('ok')

        conf.env.INCLUDES_PARALELO = includes_check
        conf.env.LIBPATH_PARALELO = libs_check
        conf.env.LIB_PARALELO = ['paralelo']

    except:
        conf.end_msg('Not found', 'RED')