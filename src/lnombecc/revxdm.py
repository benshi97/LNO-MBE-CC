from __future__ import annotations

from logging import getLogger

from ase.calculators.aims import Aims, AimsProfile

LOGGER = getLogger(__name__)


class revXDM(Aims):
    """
    FHI-aims calculator configured for revXDM calculations with enforced settings.

    This calculator pins the following parameters and does not allow them
    to be overridden:

        - species_dir : lightdense species defaults
        - xc          : "b86bpbe-50"
        - xdm         : "0.74 1.72"

    Either a command or an AimsProfile must be provided to define the
    FHI-aims executable.

    Parameters
    ----------
    command : str, optional
        Path to the FHI-aims executable. If provided, an AimsProfile
        will be constructed internally.
    profile : AimsProfile, optional
        Pre-constructed AimsProfile defining the FHI-aims execution
        environment. Exactly one of `command` or `profile` must be provided.
    species_dir_lightdense : str
        Path to the lightdense species defaults directory.
    directory : str | Path, optional
        Working directory for the calculation, by default ".".
    **kwargs
        Additional FHI-aims control parameters. The keywords
        "xc", "xdm", and "species_dir" are not allowed and will
        raise an error if provided.

    Returns
    -------
    None
    """

    def __init__(
        self,
        *args,
        command: str | None = None,
        profile: AimsProfile | None = None,
        species_dir_lightdense: str,
        directory=".",
        **kwargs,
    ):
        # Disallow user override of fixed keywords
        for forbidden in ("xc", "xdm", "species_dir"):
            if forbidden in kwargs:
                raise ValueError(f"{forbidden} is fixed in revXDM and cannot be set.")

        # Require exactly one of (command, profile)
        if (command is None) == (profile is None):
            raise ValueError("Provide exactly one of `command` or `profile`.")

        if profile is None:
            profile = AimsProfile(command=command)

        # Enforced parameters
        kwargs["species_dir"] = species_dir_lightdense
        kwargs["xc"] = "b86bpbe-50"
        kwargs["xdm"] = "0.74 1.72"

        super().__init__(*args, profile=profile, directory=directory, **kwargs)
