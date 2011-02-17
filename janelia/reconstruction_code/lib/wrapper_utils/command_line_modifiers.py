#!/usr/bin/env python
"""Provide support modifying command line args for a binary executable.

Commonly, a wrapper around a binary executable will require doing something to
the arguments of that binary, and then calling the binary with the modified
arguments as input. This module provides functions to perform common tasks of
this framework.
"""

import os
import sys

def add_negative_filter(cmd_line_args=sys.argv[1:], filter_option='-f', 
                                filter_arg_syntax='neg ARG', max_value='255'):
    """Prepend a negative filter option (invert image intensity) to a command.

    Various commands, such as watershed, expect high values for boundaries.
    In EM microscopy images, however, boundaries are typically dark. This 
    function assumes that the binary executable for a command contains an 
    option for inverting the image, and makes that option the default by 
    adding it to the command line arguments.
    """
    try:
        filter_flag_index = cmd_line_args.index('-f')
    except ValueError:
        try:
            filter_flag_index = cmd_line_args.index('-filter')
        except ValueError:
            filter_flag_index = -1
    if filter_flag_index >= 0:
        cmd_line_args[filter_flag_index+1] = \
                                filter_arg_syntax.replace('ARG', max_value) + \
                                cmd_line_args[filter_flag_index+1]
    else:
        cmd_line_args.insert(0, filter_arg_syntax.replace('ARG', max_value))
        cmd_line_args.insert(0, filter_option)
    return cmd_line_args

def run_command(binary, argv, modifier_function=add_negative_filter):
    """Run a command after first modifying the command line arguments.

    If no arguments are provided, it just runs the binary executable with no
    arguments (presumably to output usage message).
    """
    if len(argv) == 1:
        os.system(binary)
    else:
        argv = modifier_function(argv[1:])
        cmd = binary + ' ' + ' '.join(map(lambda x: '"'+x+'"', argv))
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

if __name__ == '__main__':
    run_command(sys.argv[1], sys.argv[2:], add_negative_filter)
