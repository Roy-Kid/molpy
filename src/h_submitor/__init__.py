from functools import wraps

def local(func=None):
    if func is None:
        def decorator(f):
            return f
        return decorator
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper
