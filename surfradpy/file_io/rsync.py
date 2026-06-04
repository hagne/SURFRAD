import pathlib as pl
import subprocess


def rsync_folder(p2fld_src,
                 p2fld_dst,
                 reporter=None,
                 verbose=True,
                 exclude=('.nfs*',)):
    """Sync folder contents with rsync."""
    if reporter is not None and reporter.clean <= 0:
        if verbose:
            print('no files processed => rsync skipped.')
        return None

    if verbose:
        print('starting rsync', end='...')

    p2fld_rsync_src = pl.Path(p2fld_src)
    p2fld_rsync_dst = pl.Path(p2fld_dst)
    p2fld_rsync_dst.mkdir(parents=True, exist_ok=True)

    cmd = ['rsync', '-av']
    cmd += [f'--exclude={pattern}' for pattern in exclude]
    cmd += [
        p2fld_rsync_src.as_posix().rstrip('/') + '/',
        p2fld_rsync_dst.as_posix().rstrip('/') + '/',
    ]

    try:
        out = subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        if reporter is not None:
            reporter.errors_increment(reporter.clean)
            reporter.wrapup()
        raise

    if verbose:
        print('done')

    return out
