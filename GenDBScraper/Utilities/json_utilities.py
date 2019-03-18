import json
class JSONEncoder(json.JSONEncoder):
    """ """
    """ Credits: https://stackoverflow.com/questions/33061302/dictionary-of-panda-dataframe-to-json """
    def default(self, obj):
        if hasattr(obj, 'to_json'):
            return obj.to_json()
        return json.JSONEncoder.default(self, obj)
