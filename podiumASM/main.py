#!/usr/bin/env python3
import click
import podiumASM
from podiumASM.snakeWrapper import *


@click.group(help=click.secho(podiumASM.description_tools, fg='green', nl=False), context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800},
             invoke_without_command=True, no_args_is_help=True)
@click.option('--restore', '-r', is_flag=True, required=False, default=False, show_default=True, help='Restore installation mode (need root or sudo)')
@click.version_option(podiumASM.__version__, '--version', '-v')
@click.pass_context
def main(ctx, restore):
    if ctx.invoked_subcommand is None and restore and check_privileges():
        if INSTALL_MODE.exists():
            INSTALL_MODE.unlink(missing_ok=False)
            click.secho(f"\n    Remove installation mode, now run:\n    {package_name()} install_local or install_cluster\n\n", fg="yellow")
        else:
            click.secho(f"\n    No reset need, {package_name()} not install !!!!!\n    Please run: {package_name()} install_local or install_cluster !!!!\n\n", fg="red")
    pass


# Hack for build docs with unspecified install
args = str(sys.argv)
if "sphinx" in args:
    main.add_command(podiumASM.run_cluster)
    main.add_command(podiumASM.edit_cluster_config)
    main.add_command(podiumASM.create_config)
    main.add_command(podiumASM.edit_tools)
    main.add_command(podiumASM.run_local)
    main.add_command(podiumASM.install_cluster)
    main.add_command(podiumASM.install_local)
    main.add_command(podiumASM.test_install)
else:
    mode = get_install_mode()
    if mode == "cluster":
        main.add_command(podiumASM.test_install)
        main.add_command(podiumASM.run_cluster)
        main.add_command(podiumASM.edit_cluster_config)
        main.add_command(podiumASM.create_config)
        main.add_command(podiumASM.edit_tools)
    elif mode == "local":
        main.add_command(podiumASM.test_install)
        main.add_command(podiumASM.run_local)
        main.add_command(podiumASM.create_config)
    else:
        main.add_command(podiumASM.install_cluster)
        main.add_command(podiumASM.install_local)


if __name__ == '__main__':
    main()
