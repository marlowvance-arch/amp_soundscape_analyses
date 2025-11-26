"""
naming_helpers.py
Utilities for generating standardized output filenames
based on naming rules in config.yml.
"""

from datetime import date
from typing import Iterable, Union

def _date_range_string(start: Union[str, date], end: Union[str, date]) -> str:
    if isinstance(start, date):
        start = start.strftime("%Y%m%d")
    if isinstance(end, date):
        end = end.strftime("%Y%m%d")
    return f"{start}-{end}"

def name_single_site(site: Union[int, str], test: str,
                     start_date: Union[str, date],
                     end_date: Union[str, date],
                     ext: str = "csv") -> str:
    site_str = str(site)
    dr = _date_range_string(start_date, end_date)
    fname = f"site_{site_str}_{test}_{dr}"
    return f"{fname}.{ext}"

def name_single_site_aggregated(site: Union[int, str], test: str, agg: str,
                                start_date: Union[str, date],
                                end_date: Union[str, date],
                                ext: str = "csv") -> str:
    site_str = str(site)
    dr = _date_range_string(start_date, end_date)
    fname = f"site_{site_str}_{test}_{agg}_{dr}"
    return f"{fname}.{ext}"

def name_multi_site(sites: Iterable[Union[int, str]], test: str,
                    start_date: Union[str, date],
                    end_date: Union[str, date],
                    ext: str = "csv") -> str:
    site_list = "_".join(str(s) for s in sites)
    dr = _date_range_string(start_date, end_date)
    fname = f"sites_{site_list}_{test}_{dr}"
    return f"{fname}.{ext}"

def name_all_sites(test: str,
                   start_date: Union[str, date],
                   end_date: Union[str, date],
                   ext: str = "csv") -> str:
    dr = _date_range_string(start_date, end_date)
    fname = f"all_{test}_{dr}"
    return f"{fname}.{ext}"
