class CookieCutter:
    """
    Cookie Cutter wrapper
    """

    @staticmethod
    def get_default_threads() -> int:
        return int("1")

    @staticmethod
    def get_default_mem_mb() -> int:
        return int("1024")

    @staticmethod
    def get_log_dir() -> str:
        return "cluster_logs"

    @staticmethod
    def get_default_queue() -> str:
        return ""

    @staticmethod
    def get_log_status_checks() -> bool:
        return "False" == "True"

    @staticmethod
    def get_latency_wait() -> int:
        return int("45")

    @staticmethod
    def get_max_qstat_checks() -> int:
        return int("3")

    @staticmethod
    def get_time_between_qstat_checks() -> float:
        return float("60")


