"""Event catalog for the unified Moran's I pipeline.

Defines ~24 event days spanning major storms, moderate storms, quiet days,
and substorm-active periods. Verifies data availability across all three sources.
"""

from dataclasses import dataclass, field
from pathlib import Path

import loader_tec as tec
import loader_supermag as sm
import loader_swmf as swmf


@dataclass
class Event:
    name: str
    year: int
    month: int
    day: int
    category: str  # 'major_storm', 'moderate_storm', 'quiet', 'substorm_active'
    swmf_event_id: int = None
    dst_min: float = None
    notes: str = ''
    has_swmf: bool = False
    has_tec: bool = False
    has_supermag: bool = False
    n_supermag_stations: int = 0

    @property
    def date_str(self):
        return f"{self.year:04d}-{self.month:02d}-{self.day:02d}"

    @property
    def slug(self):
        return self.name.lower().replace(' ', '_').replace("'", '')


EVENTS = [
    # === Triple overlap: SWMF + SuperMAG + TEC ===
    Event('Halloween 2003', 2003, 10, 29, 'major_storm',
          swmf_event_id=1, dst_min=-383,
          notes='Most intense storm in SWMF archive'),
    Event('Aug 2001', 2001, 8, 31, 'major_storm',
          swmf_event_id=3, dst_min=-105),
    Event('Aug 2005', 2005, 8, 31, 'moderate_storm',
          swmf_event_id=4, dst_min=-131),
    Event('Dec 2006', 2006, 12, 14, 'major_storm',
          swmf_event_id=2, dst_min=-146),

    # === Dual: SWMF + SuperMAG (no TEC on disk for 2010-2011) ===
    Event('Apr 2010', 2010, 4, 5, 'moderate_storm',
          swmf_event_id=5, dst_min=-72),
    Event('Aug 2011', 2011, 8, 5, 'moderate_storm',
          swmf_event_id=6, dst_min=-115),

    # === Dual: SuperMAG + TEC (2012-2016 processed pool) ===
    # Major storms
    Event('St Patricks 2015', 2015, 3, 17, 'major_storm', dst_min=-223),
    Event('Jun 2015', 2015, 6, 22, 'major_storm', dst_min=-204),
    Event('St Patricks 2013', 2013, 3, 17, 'moderate_storm', dst_min=-132),
    Event('Mar 2012', 2012, 3, 9, 'moderate_storm', dst_min=-133),
    Event('Oct 2012', 2012, 10, 1, 'moderate_storm', dst_min=-119),

    # Moderate storms
    Event('Feb 2014', 2014, 2, 19, 'moderate_storm', dst_min=-116),
    Event('Sep 2014', 2014, 9, 12, 'moderate_storm', dst_min=-88),
    Event('Dec 2015', 2015, 12, 20, 'moderate_storm', dst_min=-155),

    # Quiet days
    Event('Quiet Winter 2014', 2014, 1, 15, 'quiet', dst_min=-5),
    Event('Quiet Spring 2014', 2014, 3, 22, 'quiet', dst_min=-8),
    Event('Quiet Summer 2014', 2014, 7, 20, 'quiet', dst_min=-3),
    Event('Quiet Fall 2013', 2013, 9, 22, 'quiet', dst_min=-10),

    # Substorm-active but non-storm
    Event('Active Jan 2014', 2014, 1, 9, 'substorm_active', dst_min=-30),
    Event('Active Nov 2013', 2013, 11, 9, 'substorm_active', dst_min=-45),

    # === Fallback major storms from 1998-2006 (raw IONEX on disk) ===
    Event('Bastille Day 2000', 2000, 7, 15, 'major_storm', dst_min=-301),
    Event('Nov 2003', 2003, 11, 20, 'major_storm', dst_min=-422),
    Event('Nov 2004', 2004, 11, 8, 'major_storm', dst_min=-373),
    Event('Jul 2004', 2004, 7, 27, 'major_storm', dst_min=-170),
]


def verify_availability():
    """Check each event's data availability across all three sources."""
    for ev in EVENTS:
        ev.has_swmf = ev.swmf_event_id is not None and \
                       swmf.find_mag_file(ev.swmf_event_id) is not None

        ev.has_tec = tec.check_date_available(ev.year, ev.month, ev.day)

        stations = sm.available_stations_for_date(ev.year, ev.month)
        ev.n_supermag_stations = len(stations)
        ev.has_supermag = ev.n_supermag_stations >= 5


def get_verified_events(min_sources=2):
    """Return events with at least min_sources available."""
    verify_availability()
    counts = []
    for ev in EVENTS:
        n = sum([ev.has_swmf, ev.has_tec, ev.has_supermag])
        if n >= min_sources:
            counts.append((ev, n))
    return [ev for ev, _ in sorted(counts, key=lambda x: -x[1])]


def print_catalog():
    """Print formatted event catalog with availability."""
    verify_availability()
    print(f"{'Event':<25} {'Date':<12} {'Category':<18} "
          f"{'SWMF':>5} {'TEC':>5} {'SM(N)':>6} {'Sources':>8}")
    print('-' * 90)
    for ev in EVENTS:
        n_src = sum([ev.has_swmf, ev.has_tec, ev.has_supermag])
        print(f"{ev.name:<25} {ev.date_str:<12} {ev.category:<18} "
              f"{'  Y' if ev.has_swmf else '  -':>5} "
              f"{'  Y' if ev.has_tec else '  -':>5} "
              f"{ev.n_supermag_stations:>6} "
              f"{n_src:>8}")


if __name__ == '__main__':
    print_catalog()
