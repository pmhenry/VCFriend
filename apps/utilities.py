# utilities

def unquote(string): 

	if "\"" in string:
		string = string.replace('\"', '')

	if "\'" in string:
		string = string.replace('\'', '')

	return string 


#  matching algorithm, returns either '1' (match) or '0' (no match) 
def pat_match(pattern, text):  

		pattern = pattern.split(",")

		for i in range(len(pattern)):  #  compares pattern and text by individual characters
			match = True
			if pattern[i] == text[i] or pattern[i] == 'N':  #  loop continues and match remains equal to true if pattern matches  
				continue                                    #  text at given position or pattern at the given position is 'N'
			elif pattern[i] == 'Y' and text[i] != '.':
				continue                               
			else:  #  if pattern does not equal match at given position, match is set to false and breaks out of the loop
				match = False
				break		
		if match == True:
			m = 1
		else:
			m = 0

		return m
